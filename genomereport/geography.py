from pathlib import Path

GEOGRAPHY_DIR = (
    Path(__file__).parent
    / "data"
)

SPAIN_CCAA_GEOJSON = (
    GEOGRAPHY_DIR
    / "spain_ccaa.geojson"
)

AUTONOMOUS_COMMUNITY_MAP = {

    "Andalucia": "Andalucía",
    "Aragon": "Aragón",
    "Asturias": "Principado de Asturias",
    "Illes Balears": "Illes Balears",
    "Canarias": "Canarias",
    "Cantabria": "Cantabria",
    "Castilla-La Mancha": "Castilla-La Mancha",
    "Castilla y Leon": "Castilla y León",
    "Catalunya": "Cataluña/Catalunya",
    "Ceuta": "Ciudad de Ceuta",
    "Extremadura": "Extremadura",
    "Galicia": "Galicia",
    "La Rioja": "La Rioja",
    "Comunidad de Madrid": "Comunidad de Madrid",
    "Melilla": "Ciudad de Melilla",
    "Region de Murcia": "Región de Murcia",
    "Navarra": "Comunidad Foral de Navarra",
    "Pais Vasco": "País Vasco/Euskadi",
    "Comunitat Valenciana": "Comunitat Valenciana",
}

MAP_LABEL_OFFSETS = {

    "Canarias": (0.0, 0.3),
    "Illes Balears": (0.2, 0.0),
    "Ciudad de Ceuta": (0.2, -0.1),
    "Ciudad de Melilla": (0.2, -0.1),

}