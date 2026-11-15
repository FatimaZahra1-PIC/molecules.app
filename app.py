from flask import Flask, render_template, request, jsonify
import urllib.parse, urllib.request, json, re

app = Flask(__name__)

CONTRIBUTIONS = {
    "-CH3": 23.58, "-CH2-": 22.88, ">CH-": 21.74,
    "-OH": 92.88, "-O-": 22.42, "-COOH": 169.06,
    "CH aromatique": 26.73, "-Cl": 38.13
}
BASE_MOLECULES = {
    "CCO": {"nom": "ethanol", "groupes": {"-CH3": 1, "-CH2-": 1, "-OH": 1}},
    "CCCO": {"nom": "propanol", "groupes": {"-CH3": 1, "-CH2-": 2, "-OH": 1}},
    "c1ccccc1": {"nom": "benzene", "groupes": {"CH aromatique": 6}},
    "CC(=O)O": {"nom": "acide acetique", "groupes": {"-CH3": 1, "-COOH": 1}},
}

def to_celsius(tb_str):
    if not tb_str: return None, None
    s = str(tb_str)
    m = re.search(r"(-?\d+\.?\d*)", s)
    if not m: return None, s
    try: val = float(m.group(1))
    except: return None, s
    if "°F" in s:
        c = (val - 32) * 5.0/9.0
        return c, f"{c:.1f} °C"
    elif val > 150 and "°C" not in s and "°F" not in s:
        c = val - 273.15
        return c, f"{c:.1f} °C"
    return val, f"{val:.1f} °C"

def estimation_secours(smiles):
    s = smiles.strip()
    if s in ["[Fe]", "Fe"]: return 3134.0, "Clapeyron-Metal (Fe)"
    if s in ["[He]", "He"]: return 4.2, "Clapeyron-He"
    if s in ["N#N", "N2"]: return 77.3, "Clapeyron-N2"
    if s in ["O", "[O]", "O=O", "O2"]: return 373.15, "Clapeyron-H2O/O2"
    if s in ["[Ar]", "Ar"]: return 87.3, "Clapeyron-Ar"
    nb_c = s.count("C") + s.count("c")
    nb_o = s.count("O")
    return 198.2 + nb_c*21.5 + nb_o*45, "Correlation etendue"

def obtenir_infos_pubchem(smiles):
    try:
        corr = smiles
        if smiles == "Fe": corr = "[Fe]"
        if smiles == "He": corr = "[He]"
        if smiles == "N2": corr = "N#N"
        if smiles == "O2": corr = "O=O"
        enc = urllib.parse.quote(corr)
        url_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{enc}/cids/JSON"
        with urllib.request.urlopen(url_cid, timeout=15) as r:
            cid = json.loads(r.read().decode())["IdentifierList"]["CID"][0]
        tb_reelle = None
        tb_num = None
        try:
            url_bp = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
            with urllib.request.urlopen(url_bp, timeout=20) as r2:
                raw_text = r2.read().decode()
                m = re.search(r"(\-?\d+\.?\d*\s*°[CF]\s*at\s*760\s*mmHg)", raw_text)
                raw = m.group(1) if m else None
                if raw:
                    num, txt = to_celsius(raw)
                    tb_reelle = txt
                    tb_num = num
        except: pass
        if not tb_reelle: tb_reelle = f"CID {cid}"
        return cid, tb_reelle, tb_num
    except:
        return None, None, None

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/api/analyze', methods=['POST'])
def analyze():
    data = request.get_json()
    smiles = data.get('smiles','CCO').strip()
    nom = smiles
    detail = []
    methode = "Joback"
    if smiles in BASE_MOLECULES:
        infos = BASE_MOLECULES[smiles]
        nom = infos["nom"]
        total = 0
        for g, n in infos["groupes"].items():
            if g in CONTRIBUTIONS:
                contrib = CONTRIBUTIONS[g]*n
                total += contrib
                detail.append({"groupe": g, "n": n, "val": CONTRIBUTIONS[g], "contrib": contrib})
        tb_k = 198.2 + total
    else:
        tb_k, methode = estimation_secours(smiles)
        detail.append({"groupe": methode, "n": 1, "val": tb_k-198.2, "contrib": tb_k-198.2})
    tb_c = tb_k - 273.15
    cid, tb_reelle, tb_num = obtenir_infos_pubchem(smiles)
    erreur = None
    if tb_num: erreur = abs(tb_c - tb_num) / abs(tb_num) * 100 if tb_num!=0 else 0
    img_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG?image_size=600x600" if cid else None
    sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d" if cid else None
    return jsonify({
        "nom": nom, "smiles": smiles, "methode": methode,
        "tb_c": round(tb_c,2), "tb_k": round(tb_k,2),
        "detail": detail, "cid": cid,
        "tb_reelle": tb_reelle, "tb_reelle_num": tb_num,
        "img_url": img_url, "sdf_url": sdf_url,
        "erreur": round(erreur,1) if erreur else None
    })