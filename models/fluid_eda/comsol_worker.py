"""Private JSONL worker: COMSOL's Java VM is owned by this process only."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time
import traceback

import numpy as np
from scipy.integrate import cumulative_trapezoid

from .comsol_backend import COMSOL_EXECUTABLE, PROTOCOL, TEMPLATE, normalize_parameters
from models.comsol_pipe.steel_properties import thermal_conductivity_w_m_k
from models.comsol_pipe.helium_properties import saturation_pressure
from models.one_k_pot.he4_viscosity import gas_viscosity


def send(kind, request_id=None, **payload):
    print(PROTOCOL + json.dumps({"type": kind, "id": request_id, **payload},
                                ensure_ascii=False, allow_nan=False), flush=True)


def initial_temperature_profile(inputs):
    """Temperature guess only; no advective thermal boundary-value solve."""
    ta, tb = inputs["temperature_a"], inputs["temperature_b"]
    if inputs.get("wall_boundary") == "prescribed":
        knots, temperatures = wall_knots(inputs)
        x = np.unique(np.r_[np.linspace(0, 1, 4001), knots / inputs["length_m"]])
        x = x[np.r_[True, np.diff(x) > 1e-12]]
        x[-1] = 1.0
        temperature = np.interp(x * inputs["length_m"], knots, temperatures)
    elif ta == tb:
        x = np.linspace(0, 1, 501)
        temperature = np.full_like(x, ta)
    elif min(ta, tb) < 1 or max(ta, tb) > 300:
        # A linear initial guess avoids invoking the bounded SS304 fit outside
        # its source range. This does not alter the COMSOL material law.
        x = np.linspace(0, 1, 4001)
        temperature = np.linspace(ta, tb, 4001)
    else:
        temperature = np.linspace(ta, tb, 4001)
        integral = cumulative_trapezoid(thermal_conductivity_w_m_k(temperature), temperature, initial=0)
        x = integral / integral[-1]
    return x, temperature


def initial_profiles(inputs):
    """Conductive wall and pressure-squared guesses for either flow direction."""
    x, temperature = initial_temperature_profile(inputs)
    integrand = np.asarray([gas_viscosity(float(t)) * t for t in temperature])
    resistance = cumulative_trapezoid(integrand, x, initial=0)
    pa, pb = inputs["pressure_a"], inputs["pressure_b"]
    pressure = np.sqrt(np.maximum(pa**2 + (pb**2 - pa**2) * resistance / resistance[-1], min(pa, pb)**2))
    diameter = inputs["outer_diameter_m"] - 2 * inputs["wall_thickness_m"]
    mdot = np.pi * diameter**4 * (pa**2 - pb**2) / (256 * (8.314462618 / .004002602) * inputs["length_m"] * resistance[-1])
    return x, temperature, pressure, float(mdot)


def wall_knots(inputs):
    """Strictly increasing axial knots of the continuous imposed wall profile."""
    positions = [0.0]
    temperatures = [inputs["wall_profile"][0]["T_start_K"]]
    for section in inputs["wall_profile"]:
        positions.append(positions[-1] + section["length_m"])
        temperatures.append(section["T_end_K"])
    return np.asarray(positions), np.asarray(temperatures)


def wall_mesh_positions(inputs):
    """Resolve every wall-temperature slope change and its radial length scale.

    Mirrored points keep this mapped edge distribution orientation independent.
    The nominal count supplies the background grid; local refinement adds cells.
    """
    knots, _ = wall_knots(inputs)
    length = inputs["length_m"]
    positions = list(np.linspace(0, 1, inputs["axial_elements"] + 1))
    positions.extend(knots / length)
    distances = inputs["wall_thickness_m"] * np.array([.1, .3, 1, 3]) / length
    for knot in knots / length:
        positions.extend(knot + distances)
        positions.extend(knot - distances)
    positions = np.asarray(positions)
    positions = positions[(positions >= 0) & (positions <= 1)]
    positions = np.unique(np.round(np.r_[positions, 1 - positions], 13))
    return positions


def checked_energy_balance(values):
    """Reject computed states whose integrated stationary energy does not close."""
    residual = values["heat_a"] + values["heat_b"] + values["heat_environment_W"] + values["flow_energy_b"] - values["flow_energy_a"]
    dominant = max(abs(values[key]) for key in ("heat_a", "heat_b", "heat_environment_W", "flow_energy_a", "flow_energy_b"))
    tolerance = max(1e-8, 1e-3 * dominant)
    if abs(residual) > tolerance:
        raise RuntimeError(f"COMSOL pipe energy conservation failed: residual {residual:.6g} W exceeds {tolerance:.6g} W tolerance.")
    return residual, dominant


class PipeSolver:
    def __init__(self, template, executable, cores):
        os.environ["PATH"] = str(executable.parent) + os.pathsep + os.environ["PATH"]
        import mph
        self.client = mph.Client(cores=cores, version="6.1")
        if os.environ.get("FLUID_EDA_COMSOL_PROGRESS") == "1":
            from jpype import JClass
            JClass("com.comsol.model.util.ModelUtil").showProgress(True)
        self.model = self.client.load(str(template))
        self.template = template.resolve()
        self._geometry = None
        co = self.model.java.component("comp1")
        for tag, kind in (("maxgas", "Maximum"), ("mingas", "Minimum")):
            if tag not in list(co.cpl().tags()):
                op = co.cpl().create(tag, kind)
                op.selection().named("geom1_gas_dom")
        if "intouter" not in list(co.cpl().tags()):
            op = co.cpl().create("intouter", "Integration")
            op.selection().geom("geom1", 1)
            op.selection().named("outer")
        v = co.variable("diagnostics")
        for key, expression in {
            "divU": "d(u,r)+if(r>0[m],u/r,d(u,r))+d(w,z)",
            "tauzz": "2*muHe(T)*d(w,z)-2/3*muHe(T)*divU",
            "taurz": "muHe(T)*(d(u,z)+d(w,r))",
            "Fenergy": "rhoHe*w*(cpHe*T+0.5*(u^2+w^2))-u*taurz-w*tauzz",
        }.items():
            v.set(key, expression)

    def solve(self, data, save_path=None, progress=lambda message: None):
        from jpype import JInt
        start = time.monotonic()
        model, j = self.model, self.model.java
        co = j.component("comp1")
        # Always rebuild the initial guess from this exact state. A previous
        # converged solution may have different geometry, pressure, or direction.
        for parameter, key, unit in (("L", "length_m", "m"), ("Do", "outer_diameter_m", "m"),
                                     ("wall_w", "wall_thickness_m", "m"),
                                     ("P1", "pressure_a", "Pa"), ("P2", "pressure_b", "Pa"),
                                     ("T1", "temperature_a", "K"), ("T2", "temperature_b", "K")):
            j.param().set(parameter, f"{data[key]:.17g}[{unit}]")
        j.param().set("load", "1")
        ri = data["outer_diameter_m"] / 2 - data["wall_thickness_m"]
        tolerance = min(1e-6, ri / 100, data["wall_thickness_m"] / 100, data["length_m"] / 100)
        j.param().set("epssel", f"{tolerance:.17g}[m]")
        x, temperature, pressure, mdot = initial_profiles(data)
        for tag, values in (("Tguess", temperature), ("pguess", pressure)):
            j.func(tag).set("table", [[f"{xx:.17g}", f"{yy:.17g}"] for xx, yy in zip(x, values)])
            j.func(tag).set("interp", "linear" if data["wall_boundary"] == "prescribed" else "piecewisecubic")
        j.param().set("mdot_guess", f"{mdot:.17g}[kg/s]")
        spf = co.physics("spf")
        spf.feature("init1").set("p_init", "pguess(z/L)")
        spf.feature("init1").set("u_init", ["0", "0", "2*mdot_guess*RHe*Tguess(z/L)/(pguess(z/L)*pi*Ri^2)*(1-(r/Ri)^2)"])
        co.physics("ht").feature("init1").set("Tinit", "Tguess(z/L)")
        ht = co.physics("ht")
        if data["wall_boundary"] == "prescribed":
            if "edaWallT" not in list(j.func().tags()):
                wall_function = j.func().create("edaWallT", "Interpolation")
                wall_function.set("funcname", "edaWallT")
                wall_function.set("source", "table")
                wall_function.set("argunit", "m")
                wall_function.set("fununit", "K")
                wall_function.set("interp", "linear")
                wall_function.set("extrap", "const")
                wall_function.label("User-prescribed outer-wall temperature profile")
            positions, temperatures = wall_knots(data)
            j.func("edaWallT").set("table", [[f"{xx:.17g}", f"{tt:.17g}"] for xx, tt in zip(positions, temperatures)])
            if "edaWallTemp" not in list(ht.feature().tags()):
                wall_feature = ht.create("edaWallTemp", "TemperatureBoundary", 1)
                wall_feature.selection().named("outer")
                wall_feature.set("T0", "edaWallT(z)")
                wall_feature.label("Imposed external wall temperature")
            ht.feature("edaWallTemp").active(True)
        elif "edaWallTemp" in list(ht.feature().tags()):
            ht.feature("edaWallTemp").active(False)
        # COMSOL's default minimum automatic damping is 1 (a full Newton
        # step). Piecewise imposed cryogenic wall temperatures need line search.
        nonlinear = j.sol("sol1").feature("s1").feature("fc1")
        nonlinear.set("dtech", "auto")
        nonlinear.set("minstep", 1e-5)
        nonlinear.set("initstep", 0.25)
        geometry = tuple(data[key] for key in ("length_m", "outer_diameter_m", "wall_thickness_m", "axial_elements", "radial_elements")) + (data["wall_boundary"], tuple(row["length_m"] for row in data["wall_profile"]))
        j.sol("sol1").clearSolutionData()
        if geometry != self._geometry:
            progress("COMSOL: rebuilding parameterized geometry and mapped mesh…")
            co.geom("geom1").run()
            mesh = co.mesh("mesh1")
            axial = mesh.feature("map1").feature("axial")
            if data["wall_boundary"] == "prescribed":
                positions = wall_mesh_positions(data)
                axial.set("type", "explicit")
                axial.set("explicit", positions.tolist())
                self._axial_count = len(positions) - 1
            else:
                axial.set("type", "predefined")
                axial.set("elemcount", JInt(data["axial_elements"]))
                self._axial_count = data["axial_elements"]
            mesh.feature("map1").feature("radialgas").set("numelem", JInt(data["radial_elements"]))
            mesh.feature("map1").feature("radialsteel").set("numelem", JInt(max(4, data["radial_elements"] // 4)))
            mesh.run()
            self._geometry = geometry
        progress(f"COMSOL: solving pipe, {data['pressure_a']:.6g} → {data['pressure_b']:.6g} Pa, {data['temperature_a']:.5g} → {data['temperature_b']:.5g} K…")
        j.sol("sol1").runAll()
        if j.sol("sol1").isEmpty():
            raise RuntimeError("COMSOL returned without a computed pipe solution.")
        progress("COMSOL: extracting pressure, temperature, flow, and conservation checks…")
        result = self.extract(data)
        result["solve_seconds"] = time.monotonic() - start
        result["artifact"] = None
        if save_path:
            output = Path(save_path).resolve()
            if output == self.template or output == TEMPLATE.resolve():
                raise ValueError("Saving over the baseline COMSOL template is forbidden.")
            output.parent.mkdir(parents=True, exist_ok=True)
            j.label(output.name)
            j.comments("Fluid EDA parameterized pipe solution. See JSON result for exact inputs, numerical checks, and continuum applicability diagnostics.")
            model.save(str(output))
            result["artifact"] = str(output)
        return result

    def _sample(self, tag, expressions, radii, z, *, boundary=None):
        j = self.model.java
        if tag in list(j.result().numerical().tags()):
            j.result().numerical().remove(tag)
        ev = j.result().numerical().create(tag, "Interp")
        ev.set("expr", expressions)
        ev.set("coord", np.array([radii, z]).tolist())
        if boundary is not None:
            # Numerical Interp requires explicit entities, not a named selection.
            # Resolve these after each geometry rebuild instead of hard-coding IDs.
            entities = list(j.component("comp1").selection(boundary).entities(1))
            if not entities:
                raise RuntimeError(f"COMSOL boundary selection {boundary} is empty in {tag}.")
            ev.selection().geom(1)
            ev.selection().set(entities)
        raw = np.asarray(ev.getData(), dtype=float)
        expected = (len(expressions), 1, len(z))
        if raw.shape != expected:
            raise RuntimeError(f"COMSOL interpolation {tag} returned shape {raw.shape}; expected {expected}.")
        values = raw[:, 0, :]
        invalid = np.argwhere(~np.isfinite(values))
        if len(invalid):
            expression_index, point_index = invalid[0]
            raise RuntimeError(
                f"COMSOL interpolation {tag} returned a nonfinite value for "
                f"{expressions[expression_index]} at r={radii[point_index]:.17g} m, "
                f"z={z[point_index]:.17g} m "
                f"({boundary or 'domain'} selection)."
            )
        return values

    def _sample_profile(self, tag, expressions, radii, z, *, length, region="gas"):
        """Evaluate exact endpoints on their boundary, preserving point order."""
        radii, z = np.asarray(radii), np.asarray(z)
        values = np.empty((len(expressions), len(z)))
        ends = {"gas": ("inlet", "outlet"), "steel": ("steelends", "steelwarm")}[region]
        at_a, at_b = z == 0, z == length
        interior = ~(at_a | at_b)
        if interior.any():
            values[:, interior] = self._sample(tag + "_interior", expressions, radii[interior], z[interior])
        for mask, boundary in zip((at_a, at_b), ends):
            if mask.any():
                values[:, mask] = self._sample(tag + "_" + boundary, expressions,
                                               radii[mask], z[mask], boundary=boundary)
        return values

    def extract(self, data):
        def scalar(expression):
            value = float(self.model.evaluate(expression))
            if not np.isfinite(value):
                raise RuntimeError(f"Non-finite COMSOL result: {expression}")
            return value

        expressions = {
            "mass_in": "comp1.intin(2*pi*r*comp1.rhoHe*w)",
            "mass_out": "comp1.intout(2*pi*r*comp1.rhoHe*w)",
            "heat_a": "comp1.intcold(2*pi*r*comp1.ht.ndflux)",
            "heat_b": "comp1.intwarm(2*pi*r*comp1.ht.ndflux)",
            "heat_environment_W": "comp1.intouter(2*pi*r*comp1.ht.ndflux)",
            "flow_energy_a": "comp1.intin(2*pi*r*comp1.Fenergy)",
            "flow_energy_b": "comp1.intout(2*pi*r*comp1.Fenergy)",
            "pressure_a_mean": "comp1.intin(2*pi*r*p)/comp1.intin(2*pi*r)",
            "pressure_b_mean": "comp1.intout(2*pi*r*p)/comp1.intout(2*pi*r)",
            "max_kn": "comp1.maxgas(comp1.Kn)",
            "max_mach": "comp1.maxgas(comp1.Mach)",
            "max_re": "comp1.maxgas(comp1.Re_local)",
            "min_pressure": "comp1.mingas(p)",
            "min_temperature": "comp1.mingas(T)",
            "max_temperature": "comp1.maxgas(T)",
            "inlet_area": "comp1.intin(2*pi*r)",
        }
        values = {name: scalar(expression) for name, expression in expressions.items()}
        mass_error = abs(values["mass_in"] - values["mass_out"]) / max(abs(values["mass_in"]), abs(values["mass_out"]), 1e-18)
        energy_error, dominant_energy_flux = checked_energy_balance(values)
        if values["min_pressure"] <= 0:
            raise RuntimeError(f"Solved pipe minimum absolute pressure must be positive: {values['min_pressure']:.9g} Pa.")
        if values["min_temperature"] <= 0:
            raise RuntimeError(f"Solved pipe minimum absolute temperature must be positive: {values['min_temperature']:.9g} K.")
        if abs(values["mass_in"] - values["mass_out"]) > max(1e-15, 1e-3 * max(abs(values["mass_in"]), abs(values["mass_out"]))):
            raise RuntimeError(f"COMSOL pipe mass conservation failed: relative mismatch {mass_error:.3g}.")
        length = data["length_m"]
        ri = data["outer_diameter_m"] / 2 - data["wall_thickness_m"]
        z = np.linspace(0, length, data["profile_points"])
        # Gauss quadrature in q=(r/R)^2 integrates an area mean directly.
        nodes, weights = np.polynomial.legendre.leggauss(32)
        radius = ri * np.sqrt((nodes + 1) / 2)
        radii = np.tile(radius, len(z))
        zz = np.repeat(z, len(radius))
        sampled = self._sample_profile("eda_bulk", ["p", "T", "w", "rhoHe*w", "rhoHe*w*T"], radii, zz, length=length).reshape(5, len(z), len(radius))
        averaged = np.einsum("ijr,r->ij", sampled, weights / 2)
        axis = self._sample_profile("eda_axis", ["p", "T", "w", "Kn"], np.zeros_like(z), z, length=length)
        wall = self._sample_profile("eda_wall", ["T"], np.full_like(z, ri + data["wall_thickness_m"] / 2), z, length=length, region="steel")[0]
        wall_outer = self._sample_profile("eda_wall_outer", ["T"], np.full_like(z, data["outer_diameter_m"] / 2), z, length=length, region="steel")[0]
        profile = []
        for i, position in enumerate(z):
            profile.append({"x": float(position), "pressure": float(averaged[0, i]),
                            "temperature": float(averaged[1, i]), "velocity": float(averaged[2, i]),
                            "wall_temperature": float(wall[i]),
                            "outer_wall_temperature": float(wall_outer[i]),
                            "bulk_temperature": float(averaged[4, i] / averaged[3, i]) if abs(averaged[3, i]) > 1e-18 else float(averaged[1, i]),
                            "axis_pressure": float(axis[0, i]), "axis_temperature": float(axis[1, i]),
                            "axis_velocity": float(axis[2, i]), "knudsen": float(axis[3, i])})
            for alias, key in (("x_m", "x"), ("pressure_Pa", "pressure"), ("temperature_K", "temperature"),
                               ("wall_temperature_K", "wall_temperature"), ("velocity_m_s", "velocity"), ("Kn", "knudsen")):
                profile[-1][alias] = profile[-1][key]
        diagnostics = {**values, "relative_mass_error": mass_error,
                       "energy_residual_W": energy_error,
                       "relative_energy_error": abs(energy_error) / dominant_energy_flux if dominant_energy_flux else 0.0,
                       "continuum_baseline": True,
                       "continuum_applicable": values["max_kn"] < 0.01,
                       "max_Kn": values["max_kn"], "max_Mach": values["max_mach"], "max_Re": values["max_re"],
                       "mesh": {"axial": self._axial_count, "axial_nominal": data["axial_elements"], "radial_gas": data["radial_elements"], "radial_steel": max(4, data["radial_elements"] // 4)}}
        warnings = []
        phase_ratios = [row["pressure"] / saturation_pressure(row["temperature"])
                        for row in profile if 0.8 <= row["temperature"] <= 4.3]
        diagnostics["max_pressure_over_saturation"] = float(max(phase_ratios)) if phase_ratios else None
        diagnostics["condensation_warning"] = bool(phase_ratios and max(phase_ratios) > 1.001)
        if diagnostics["condensation_warning"]:
            warnings.append("Sampled helium pressure exceeds saturation: the single-phase ideal-gas pipe model is not physically valid there.")
        if values["max_kn"] >= 0.01:
            warnings.append("Kn ≥ 0.01: the no-slip continuum pipe model is outside its strict continuum regime.")
        if values["max_mach"] >= 0.3:
            warnings.append("Mach ≥ 0.3: the template's low-Mach flow assumption is exceeded.")
        if values["max_re"] >= 2000:
            warnings.append("Local diameter-based Reynolds number ≥ 2000: check the laminar assumption.")
        return {"backend": "COMSOL 6.1", "mass_flow": (values["mass_in"] + values["mass_out"]) / 2,
                "mass_flow_kg_s": (values["mass_in"] + values["mass_out"]) / 2,
                "pressure_in": values["pressure_a_mean"], "pressure_out": values["pressure_b_mean"],
                "temperature_in": data["temperature_a"], "temperature_out": data["temperature_b"],
                "heat_a": values["heat_a"], "heat_b": values["heat_b"],
                "heat_a_W": values["heat_a"], "heat_b_W": values["heat_b"],
                "heat_environment_W": values["heat_environment_W"],
                "profile": profile, "diagnostics": diagnostics, "warnings": warnings}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path, default=TEMPLATE)
    parser.add_argument("--executable", type=Path, default=COMSOL_EXECUTABLE)
    parser.add_argument("--cores", type=int, default=2)
    args = parser.parse_args()
    try:
        solver = PipeSolver(args.template.resolve(), args.executable.resolve(), args.cores)
        send("ready")
    except Exception as exc:
        send("error", error=f"COMSOL startup failed: {exc}", traceback=traceback.format_exc())
        return 1
    for line in sys.stdin:
        request_id = None
        try:
            request = json.loads(line)
            request_id = request["id"]
            inputs = request["inputs"]
            inputs = normalize_parameters(inputs, inputs["pressure_a"], inputs["pressure_b"], inputs["temperature_a"], inputs["temperature_b"])
            result = solver.solve(inputs, request.get("save_path"), lambda message: send("progress", request_id, message=message))
            send("result", request_id, result=result)
        except Exception as exc:
            send("error", request_id, error=f"COMSOL pipe evaluation failed: {exc}", traceback=traceback.format_exc())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
