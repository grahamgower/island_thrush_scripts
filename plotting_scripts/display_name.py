
def italicise_display_name(name):
    it_name = ["\\emph{"]
    state = 1
    for field in name.split(" "):
        if state == 1 and (field.startswith("(") or field in ("ssp", "AMNH", "ZMUC")):
                it_name.append("}")
                state = 0
        it_name.append(field)
    if state == 1:
        it_name.append("}")
    disp_name = " ".join(it_name)
    disp_name = disp_name.replace("{ ", "{")
    disp_name = disp_name.replace(" }", "}")
    return disp_name
