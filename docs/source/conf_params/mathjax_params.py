mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

macros = dict()

# general coordinates
macros["gx"] = "{\\xi}"
macros["gy"] = "{\\eta}"
macros["gz"] = "{\\zeta}"

# velocity
macros["ux"] = "{u_x}"
macros["uy"] = "{u_y}"
macros["uz"] = "{u_z}"

# cell centers
macros["pimm"] = "i-1           "
macros["pim" ] = "i-\\frac{1}{2}"
macros["pic" ] = "i             "
macros["pip" ] = "i+\\frac{1}{2}"
macros["pipp"] = "i+1           "
macros["pjmm"] = "j-1           "
macros["pjm" ] = "j-\\frac{1}{2}"
macros["pjc" ] = "j             "
macros["pjp" ] = "j+\\frac{1}{2}"
macros["pjpp"] = "j+1           "
macros["pkmm"] = "k-1           "
macros["pkm" ] = "k-\\frac{1}{2}"
macros["pkc" ] = "k             "
macros["pkp" ] = "k+\\frac{1}{2}"
macros["pkpp"] = "k+1           "

# x cell faces
macros["ximm"] = "i-\\frac{1}{2}"
macros["xim" ] = "i             "
macros["xic" ] = "i+\\frac{1}{2}"
macros["xip" ] = "i+1           "
macros["xipp"] = "i+\\frac{3}{2}"
macros["xjmm"] = "j-1           "
macros["xjm" ] = "j-\\frac{1}{2}"
macros["xjc" ] = "j             "
macros["xjp" ] = "j+\\frac{1}{2}"
macros["xjpp"] = "j+1           "
macros["xkmm"] = "k-1           "
macros["xkm" ] = "k-\\frac{1}{2}"
macros["xkc" ] = "k             "
macros["xkp" ] = "k+\\frac{1}{2}"
macros["xkpp"] = "k+1           "

# y cell faces
macros["yimm"] = "i-1           "
macros["yim" ] = "i-\\frac{1}{2}"
macros["yic" ] = "i             "
macros["yip" ] = "i+\\frac{1}{2}"
macros["yipp"] = "i+1           "
macros["yjmm"] = "j-\\frac{1}{2}"
macros["yjm" ] = "j             "
macros["yjc" ] = "j+\\frac{1}{2}"
macros["yjp" ] = "j+1           "
macros["yjpp"] = "j+\\frac{3}{2}"
macros["ykmm"] = "k-1           "
macros["ykm" ] = "k-\\frac{1}{2}"
macros["ykc" ] = "k             "
macros["ykp" ] = "k+\\frac{1}{2}"
macros["ykpp"] = "k+1           "

# z cell faces
macros["zimm"] = "i-1           "
macros["zim" ] = "i-\\frac{1}{2}"
macros["zic" ] = "i             "
macros["zip" ] = "i+\\frac{1}{2}"
macros["zipp"] = "i+1           "
macros["zjmm"] = "j-1           "
macros["zjm" ] = "j-\\frac{1}{2}"
macros["zjc" ] = "j             "
macros["zjp" ] = "j+\\frac{1}{2}"
macros["zjpp"] = "j+1           "
macros["zkmm"] = "k-\\frac{1}{2}"
macros["zkm" ] = "k             "
macros["zkc" ] = "k+\\frac{1}{2}"
macros["zkp" ] = "k+1           "
macros["zkpp"] = "k+\\frac{3}{2}"

# derivatives, continuous and discrete
macros["der"]  = ["{\\frac{\\partial #1}{\\partial #2}}", 2]
macros["dder"] = ["{\\frac{\\delta #1}{\\delta #2}}", 2]

# mesh skewness tensor
macros["mst"] = ["{\\gamma^{#1 #2}}", 2]

# interpolation and differentiation
macros["intrp"] = ["{\\overline{#1}^{#2}}", 2]
macros["diffe"] = ["{\\delta_{#2} {#1}}", 2]

# value at
macros["vat"] = ["{\\left. {#1} \\right|_{#2}}", 2]

# average
macros["ave"] = ["{\\left\\langle {#1} \\right\\rangle_{#2}}", 2]

# interpolations, arithmetic, volume, unknown
macros["dintrpa"] = ["{\\overline {#1}^{\\left( A,#2 \\right)}}", 2]
macros["dintrpv"] = ["{\\overline {#1}^{\\left( V,#2 \\right)}}", 2]
macros["dintrpu"] = ["{\\overline {#1}^{\\left( U,#2 \\right)}}", 2]

mathjax3_config = {"TeX": {"Macros": macros}}
