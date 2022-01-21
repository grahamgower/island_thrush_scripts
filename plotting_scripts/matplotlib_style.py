import matplotlib
import matplotlib.pyplot as plt

#matplotlib.use('pgf')

# set font
plt.rcParams.update(
    {
        "text.usetex": True,
        "text.latex.preamble": 
            r"""
            \usepackage[utf8]{inputenc}
            \usepackage{newunicodechar,graphicx}
            \DeclareRobustCommand{\okina}{%
              \raisebox{\dimexpr\fontcharht\font`A-\height}{%
                \scalebox{0.8}{`}%
              }%
            }
            \newunicodechar{ʻ}{\okina}
            """,
        #"pgf.preamble": r"\usepackage[utf8]{inputenc}",
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
    }
)
