"""Make the pandoc reference.docx used for the BMC manuscript: Times New Roman 12 pt, double line spacing,
2.5 cm margins, A4. Styles only; content comes from the markdown."""
import sys
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_LINE_SPACING
src, dst = sys.argv[1], sys.argv[2]
d = Document(src)
for s in d.styles:
    try:
        f = s.font
    except Exception:
        continue
    if s.type != 1 and s.type != 2:          # paragraph (1) and character (2) styles only
        continue
    f.name = "Times New Roman"
    rpr = s.element.get_or_add_rPr(); rf = rpr.find('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}rFonts')
    if rf is not None:
        for att in ["ascii", "hAnsi", "cs", "eastAsia"]:
            rf.set('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}' + att, "Times New Roman")
        for att in ["asciiTheme", "hAnsiTheme", "cstheme", "eastAsiaTheme"]:
            k = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}' + att
            if k in rf.attrib: del rf.attrib[k]
    if s.type == 1:
        pf = s.paragraph_format
        pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
        pf.space_after = Pt(0); pf.space_before = Pt(0)
    try:
        f.color.rgb = RGBColor(0, 0, 0)
    except Exception:
        pass
for name, size, bold in [("Normal", 12, False), ("Body Text", 12, False), ("First Paragraph", 12, False),
                         ("Compact", 12, False), ("Title", 14, True), ("Heading 1", 14, True),
                         ("Heading 2", 12, True), ("Heading 3", 12, True), ("Bibliography", 12, False),
                         ("Abstract", 12, False), ("Author", 12, False), ("Date", 12, False),
                         ("Image Caption", 12, False), ("Table Caption", 12, False)]:
    try:
        st = d.styles[name]
    except KeyError:
        continue
    st.font.size = Pt(size); st.font.bold = bold; st.font.italic = False if name != "Heading 3" else True
    if name.startswith("Heading") or name == "Title":
        st.paragraph_format.space_before = Pt(12); st.paragraph_format.keep_with_next = True
    if name == "Title":
        st.paragraph_format.alignment = 0
sec = d.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
for side in ["left_margin", "right_margin", "top_margin", "bottom_margin"]:
    setattr(sec, side, Cm(2.5))
d.save(dst)
print("wrote", dst)
