"""Carry component identity and solver input context through the job API."""
from .network import NetworkError


def pipe_error(error, component, pa, pb, ta, tb, *, validation=False):
    """Describe the actual attempted pipe state, without changing its values."""
    context = {"pressure_a": pa, "pressure_b": pb,
               "temperature_a": ta, "temperature_b": tb}
    reason = str(error)
    message = "COMSOL could not solve this pipe."
    hint = "Check this pipe's geometry and boundary conditions. Expand technical details for the solver message."
    if validation:
        message = reason
        hint = "Review this pipe's parameters and thermal boundaries."
        for field, value in tuple(context.items()):
            if reason.startswith(field + " "):
                quantity, port = field.split("_")
                unit = "Pa" if quantity == "pressure" else "K"
                context.update(field=field, port=port, value=value, unit=unit)
                message = reason.replace(field, f"Port {port.upper()} {quantity} (calculated {value:.6g} {unit})", 1)
                hint = ("Absolute pressure must be finite and greater than zero. "
                        "Check the connected pressure sources and pump parameters."
                        if quantity == "pressure" else
                        "Check the pipe end temperatures and connected temperature boundaries against the supported material range.")
                break
    return NetworkError(message, [component["id"]],
                        "comsol_input" if validation else "comsol_pipe_failure",
                        stage="COMSOL input validation" if validation else "COMSOL pipe solve",
                        hint=hint, context=context, details=reason)


def describe_error(error, project, *, stage="Circuit validation", code="simulation_error"):
    """Resolve stable IDs to names in the submitted schematic, including multi-device errors."""
    detail = error.as_dict() if isinstance(error, NetworkError) else {
        "message": str(error), "code": code, "component_ids": []}
    raw = project.get("components", []) if isinstance(project, dict) else []
    components = {c.get("id"): c for c in raw if isinstance(c, dict) and isinstance(c.get("id"), str)} if isinstance(raw, list) else {}
    detail["components"] = [
        {"id": cid, "label": str(components[cid].get("label") or cid),
         "type": components[cid].get("type", "")}
        for cid in detail["component_ids"] if cid in components]
    detail.setdefault("stage", stage)
    detail.setdefault("details", detail["message"])
    # Existing topology errors often start with an internal ID. Show the user's name.
    for component in detail["components"]:
        prefix = component["id"] + ":"
        if detail["message"].startswith(prefix):
            detail["message"] = component["label"] + ":" + detail["message"][len(prefix):]
            break
    detail.setdefault("hint", "Select a component below to inspect its connections and parameters."
                      if detail["components"] else "Review the circuit connections and boundary conditions. Expand technical details for more information.")
    return detail
