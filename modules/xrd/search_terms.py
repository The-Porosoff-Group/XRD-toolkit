"""Chemical names and searchable metadata for phase identification."""
import re


ELEMENT_NAMES = dict(pair.split(':') for pair in (
    'H:hydrogen He:helium Li:lithium Be:beryllium B:boron C:carbon N:nitrogen '
    'O:oxygen F:fluorine Ne:neon Na:sodium Mg:magnesium Al:aluminium Si:silicon '
    'P:phosphorus S:sulfur Cl:chlorine Ar:argon K:potassium Ca:calcium Sc:scandium '
    'Ti:titanium V:vanadium Cr:chromium Mn:manganese Fe:iron Co:cobalt Ni:nickel '
    'Cu:copper Zn:zinc Ga:gallium Ge:germanium As:arsenic Se:selenium Br:bromine '
    'Kr:krypton Rb:rubidium Sr:strontium Y:yttrium Zr:zirconium Nb:niobium '
    'Mo:molybdenum Tc:technetium Ru:ruthenium Rh:rhodium Pd:palladium Ag:silver '
    'Cd:cadmium In:indium Sn:tin Sb:antimony Te:tellurium I:iodine Xe:xenon '
    'Cs:caesium Ba:barium La:lanthanum Ce:cerium Pr:praseodymium Nd:neodymium '
    'Pm:promethium Sm:samarium Eu:europium Gd:gadolinium Tb:terbium Dy:dysprosium '
    'Ho:holmium Er:erbium Tm:thulium Yb:ytterbium Lu:lutetium Hf:hafnium '
    'Ta:tantalum W:tungsten Re:rhenium Os:osmium Ir:iridium Pt:platinum Au:gold '
    'Hg:mercury Tl:thallium Pb:lead Bi:bismuth Po:polonium At:astatine Rn:radon '
    'Fr:francium Ra:radium Ac:actinium Th:thorium Pa:protactinium U:uranium '
    'Np:neptunium Pu:plutonium Am:americium Cm:curium Bk:berkelium Cf:californium '
    'Es:einsteinium Fm:fermium Md:mendelevium No:nobelium Lr:lawrencium '
    'Rf:rutherfordium Db:dubnium Sg:seaborgium Bh:bohrium Hs:hassium Mt:meitnerium '
    'Ds:darmstadtium Rg:roentgenium Cn:copernicium Nh:nihonium Fl:flerovium '
    'Mc:moscovium Lv:livermorium Ts:tennessine Og:oganesson'
).split())
NAME_TO_ELEMENT = {name: symbol for symbol, name in ELEMENT_NAMES.items()}
NAME_TO_ELEMENT.update(aluminum='Al', sulphur='S', cesium='Cs', wolfram='W')
COMPOUND_TERMS = {
    'carbide': 'C', 'nitride': 'N', 'oxide': 'O', 'silicide': 'Si',
    'sulfide': 'S', 'sulphide': 'S', 'phosphide': 'P', 'boride': 'B',
    'hydride': 'H', 'fluoride': 'F', 'chloride': 'Cl', 'bromide': 'Br',
    'iodide': 'I', 'selenide': 'Se', 'telluride': 'Te',
}


def candidate_search_text(entry):
    """Index returned descriptions and full element names alongside formulas."""
    fields = ('formula', 'name', 'description', 'mineral', 'spacegroup',
              'system', 'mp_id', 'mp_api_id', 'cod_id', 'authors', 'journal',
              'year', 'source', 'stability')
    text = [str(entry.get(field) or '') for field in fields]
    symbols = set(re.findall(r'[A-Z][a-z]?', entry.get('formula') or ''))
    text.extend(name for name, symbol in NAME_TO_ELEMENT.items() if symbol in symbols)
    return ' '.join(text).lower()
