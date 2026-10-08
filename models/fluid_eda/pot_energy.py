"""Continuous He-4 pot enthalpy balance, excluding component-end conduction.

Pot saturation h_l and latent heat: existing Donnelly & Barenghi (1998) fits.
Normal-fluid feed h(T,P): CoolProp HEOS helium, with one constant reference
shift anchored to DB98 saturated liquid at 4.2 K. This joins two property sources;
reference alignment is not an independent validation of their differences.
Below 2.1768 K, do not extrapolate the normal-fluid EOS into superfluid helium.
Use the existing bounded DB98 liquid approximation or dilute-vapor approximation.
Hydraulic pipe/JT laws are unchanged; this evaluates refrigeration at their flow.
"""
from functools import lru_cache
import math
from models.continuous_he4 import properties as he4

REFERENCE_T_K = 4.2
NORMAL_FLUID_MIN_T_K = 2.1768
MODEL = "DB98_pot_CoolProp_feed"


@lru_cache(maxsize=1)
def reference_shift():
    from CoolProp.CoolProp import PropsSI
    return float(he4.saturated_liquid_enthalpy(REFERENCE_T_K)) - PropsSI(
        "Hmass", "T", REFERENCE_T_K, "Q", 0, "Helium")


@lru_cache(maxsize=512)
def feed_state(temperature, pressure):
    """Return h [J/kg], phase, and method with a common DB98 enthalpy zero."""
    if not all(math.isfinite(x) and x > 0 for x in (temperature, pressure)):
        raise ValueError("Feed temperature and pressure must be finite and positive.")
    if temperature < NORMAL_FLUID_MIN_T_K:
        ps = float(he4.saturation_pressure(temperature))
        if math.isclose(pressure, ps, rel_tol=1e-5):
            raise ValueError("The pre-JT feed lies on saturation; temperature and pressure alone do not determine its liquid/vapor fraction.")
        if pressure > ps:
            return float(he4.liquid_enthalpy(temperature, pressure)), "liquid", "DB98 compressed-liquid approximation"
        return float(he4.saturated_vapor_enthalpy(temperature)), "gas", "DB98 dilute-vapor approximation h(T,P)=h_sat_vapor(T)"
    from CoolProp.CoolProp import PropsSI, PhaseSI
    # Ranges follow the normal-fluid source and implementation, not the old
    # blanket 1--300 K adapter check. No extrapolation into He II is attempted.
    if temperature > 1500 or pressure > 1e9:
        raise ValueError("Normal-helium feed enthalpy supports T <= 1500 K and P <= 1e9 Pa.")
    h = PropsSI("Hmass", "T", temperature, "P", pressure, "Helium") + reference_shift()
    phase = PhaseSI("T", temperature, "P", pressure, "Helium")
    if not math.isfinite(h):
        raise ValueError("The helium feed enthalpy is not finite.")
    return float(h), phase, "CoolProp HEOS, DB98-aligned enthalpy reference"


def flash_at_pot(h, temperature, pressure):
    """Flash classification at pot pressure; never clamp quality into [0,1]."""
    hl = float(he4.saturated_liquid_enthalpy(temperature))
    hg = float(he4.saturated_vapor_enthalpy(temperature))
    fraction = (h-hl)/(hg-hl)
    if 0 <= fraction <= 1:
        return {"phase": "two_phase", "vapor_quality": fraction, "temperature_K": temperature}
    if fraction < 0:
        return {"phase": "subcooled_liquid", "vapor_quality": None, "temperature_K": None}
    tout = None
    # An EOS flash is only attempted if the gas state lies above its lower
    # temperature boundary. Low-temperature gas remains explicitly unresolved.
    try:
        from CoolProp.CoolProp import PropsSI
        raw_h = h-reference_shift()
        low_h = PropsSI("Hmass", "T", NORMAL_FLUID_MIN_T_K+1e-6, "P", pressure, "Helium")
        if raw_h >= low_h:
            candidate = float(PropsSI("T", "Hmass", raw_h, "P", pressure, "Helium"))
            if NORMAL_FLUID_MIN_T_K <= candidate <= 1500:
                tout = candidate
    except (ValueError, ImportError):
        pass
    return {"phase": "superheated_vapor", "vapor_quality": None, "temperature_K": tout}


def pot_balance(temperature, incoming, mass_out, heat_load=0.):
    """Positive net power can remove load; positive residual warms the pot.

    Q = mout*(hg-hl) - sum(min*(hin-hl)). Referencing both streams to the pot
    liquid keeps the result reference-independent even with a small mass residual.
    At steady inventory this is exactly mout*hg - sum(min*hin). Nonzero inventory
    rate is reported, not silently supplied by the fixed-pressure source.
    """
    if not math.isfinite(mass_out) or mass_out < 0 or not math.isfinite(heat_load) or heat_load < 0:
        raise ValueError("Pot outflow and applied heat load must be finite and nonnegative.")
    for stream in incoming:
        if not math.isfinite(stream["mass_flow_kg_s"]) or stream["mass_flow_kg_s"] < 0 or not math.isfinite(stream["enthalpy_J_kg"]):
            raise ValueError("Incoming pot streams require finite enthalpy and nonnegative mass flow.")
    hl = float(he4.saturated_liquid_enthalpy(temperature))
    hg = float(he4.saturated_vapor_enthalpy(temperature))
    mass_in = sum(s["mass_flow_kg_s"] for s in incoming)
    incoming_load = sum(s["mass_flow_kg_s"]*(s["enthalpy_J_kg"]-hl) for s in incoming)
    evaporation = mass_out*(hg-hl)
    net = evaporation-incoming_load
    return {"available": True, "model": MODEL,
            "liquid_enthalpy_J_kg": hl, "vapor_enthalpy_J_kg": hg,
            "latent_heat_J_kg": hg-hl, "mass_in_kg_s": mass_in,
            "mass_out_kg_s": mass_out, "inventory_rate_kg_s": mass_in-mass_out,
            "relative_mass_residual": abs(mass_in-mass_out)/max(mass_in+mass_out,1e-16),
            "evaporation_power_W": evaporation, "incoming_enthalpy_load_W": incoming_load,
            "net_cooling_power_W": net, "applied_heat_load_W": heat_load,
            "energy_residual_W": heat_load-net, "incoming_streams": incoming,
            "circulating": mass_in > 0 and mass_out > 0,
            "solid_conduction_included": False,
            "enthalpy_reference": "DB98 liquid ground-state reference; CoolProp aligned at saturated liquid 4.2 K"}


def apply_pot_balances(result, project):
    """Use accepted branch flows; never invoke an imaginary gas-thermalizing bath."""
    nodes = {n["id"]:n for n in result["nodes"]}
    components = {c["id"]:c for c in project["components"]}
    for component in project["components"]:
        params = component.get("params", {})
        if component["type"] != "source" or params.get("mode") != "saturation":
            continue
        cid = component["id"]
        row = result["components"][cid]
        node = nodes[row["node"]]
        temperature = node["temperature"]
        incoming, mass_out, changes = [], 0., []
        notes = []
        try:
            for edge in result["edges"]:
                nid = node["id"]
                if nid not in (edge["node_a"],edge["node_b"]):
                    continue
                outward = edge["mass_flow"] * (1 if edge["node_a"] == nid else -1)
                if outward >= 0:
                    mass_out += outward
                    continue
                upstream = nodes[edge["node_a"] if edge["node_b"] == nid else edge["node_b"]]
                if edge["type"] == "jt":
                    if any(components.get(sid,{}).get("params",{}).get("mode")=="saturation" for sid in upstream.get("source_ids",[])):
                        h = float(he4.saturated_vapor_enthalpy(upstream["temperature"]))
                        phase, method = "gas", "DB98 saturated-vapor source"
                    else:
                        h, phase, method = feed_state(upstream["temperature"],upstream["pressure"])
                    flash = flash_at_pot(h,temperature,node["pressure"])
                    changes.append((edge,h,flash))
                    if phase not in ("gas", "supercritical_gas"):
                        notes.append("The feed is liquid/dense fluid; current pipe and restriction mass flow still uses gas hydraulics. Cooling power at this flow is conditional, not a validated liquid/two-phase circuit solution.")
                    if "approximation" in method:
                        notes.append(method+" is used for the low-temperature feed enthalpy.")
                else:
                    # Pipe/pump models deliver gas at the specified receiving
                    # port state. At a saturated pot boundary that is vapor.
                    h = float(he4.saturated_vapor_enthalpy(temperature))
                    phase, method = "saturated_vapor", "DB98; gas at prescribed pot-side port state"
                    flash = {"phase":"saturated_vapor","vapor_quality":1.,"temperature_K":temperature}
                incoming.append({"component_id":edge["id"],"mass_flow_kg_s":-outward,
                                 "enthalpy_J_kg":h,"upstream_phase":phase,"method":method,
                                 "state_at_pot_pressure":flash})
            balance = pot_balance(temperature,incoming,mass_out,float(params.get("heat_load_W",0.)))
            if not balance["circulating"]:
                notes.append("No continuous inlet/outlet circulation. A zero heat imbalance alone does not establish continuous-pot equilibrium.")
            balance["notes"] = list(dict.fromkeys(notes))
            import CoolProp
            balance["property_sources"] = {"pot":"Donnelly & Barenghi (1998), Tables 7.5, 16.2, 17.2",
                "feed":"CoolProp HEOS helium", "CoolProp_version":CoolProp.__version__,
                "reference_temperature_K":REFERENCE_T_K,"reference_shift_J_kg":reference_shift()}
            # Only update the plotted JT exit state after all property calls
            # succeeded; mass-flow constitutive laws and solved nodes are intact.
            for edge,h,flash in changes:
                edge["hydraulic_temperature_model"] = "ideal_gas_isothermal_resistance"
                edge["thermodynamic_model"] = "real_helium_isenthalpic_pot_feed"
                edge["enthalpy_in_J_kg"] = edge["enthalpy_out_J_kg"] = h
                edge["outlet_fluid_temperature_K"] = flash["temperature_K"]
                if edge["mass_flow"] >= 0:
                    edge["temperature_out"] = flash["temperature_K"]
                else:
                    edge["temperature_in"] = flash["temperature_K"]
                    edge["temperature_out"] = nodes[edge["node_b"]]["temperature"]
                edge["flash_at_pot"] = flash
                edge["warnings"] = [w for w in edge.get("warnings",[]) if not w.startswith("Ideal-gas isenthalpic")]
                edge["warnings"].append("Outlet thermodynamics uses helium enthalpy; the restriction flow law remains a gas approximation.")
                prefix = f"{edge.get('label',edge['id'])}: Ideal-gas isenthalpic"
                result["warnings"] = [w for w in result["warnings"] if not w.startswith(prefix)]
                result["warnings"].append(f"{edge.get('label',edge['id'])}: "+edge["warnings"][-1])
                if edge.get("profile"):
                    destination = -1 if edge["mass_flow"] >= 0 else 0
                    edge["profile"][destination]["temperature"] = flash["temperature_K"]
                    for point in edge["profile"]:
                        point["wall_temperature"] = None
            row["pot_energy"] = balance
        except (ValueError, ArithmeticError, ImportError) as error:
            row["pot_energy"] = {"available":False,"model":MODEL,"message":str(error),"solid_conduction_included":False}
            notes.append("Pot energy balance unavailable: "+str(error))
        for note in dict.fromkeys(notes):
            result["warnings"].append(f"{row.get('label',cid)}: {note}")
    return result
