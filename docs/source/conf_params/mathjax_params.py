mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
mathjax3_config = {
    "TeX": {
        "Macros": {
            "der": ["{\\frac{\\partial #1}{\\partial #2}}", 2], # derivative
            "dder": ["{\\frac{\\delta #1}{\\delta #2}}", 2],    # discrete derivative
            "mst": ["{\\gamma^{#1 #2}}", 2], # mesh skewness tensor
            "gx": ["{\\xi}"],
            "gy": ["{\\eta}"],
            "gz": ["{\\zeta}"],
            "ux": ["{u_x}"],
            "uy": ["{u_y}"],
            "uz": ["{u_z}"],
            "intrp": ["{\\overline{#1}^{#2}}", 2], # interpolation
            "diffe": ["{\\delta_{#2} {#1}}", 2],   # differentiation
            "vat": ["{\\left. {#1} \\right|_{#2}}", 2], # value at
            "ave": ["{\\left\\langle {#1} \\right\\rangle_{#2}}", 2],
            ## indices, pressure, x-face, y-face in two directions
            # p-i
            "pimm": ["i-1           "],
            "pim":  ["i-\\frac{1}{2}"],
            "pic":  ["i             "],
            "pip":  ["i+\\frac{1}{2}"],
            "pipp": ["i+1           "],
            # p-j
            "pjmm": ["j-1           "],
            "pjm":  ["j-\\frac{1}{2}"],
            "pjc":  ["j             "],
            "pjp":  ["j+\\frac{1}{2}"],
            "pjpp": ["j+1           "],
            # p-k
            "pkmm": ["k-1           "],
            "pkm":  ["k-\\frac{1}{2}"],
            "pkc":  ["k             "],
            "pkp":  ["k+\\frac{1}{2}"],
            "pkpp": ["k+1           "],
            # x-i
            "ximm": ["i-\\frac{1}{2}"],
            "xim":  ["i             "],
            "xic":  ["i+\\frac{1}{2}"],
            "xip":  ["i+1           "],
            "xipp": ["i+\\frac{3}{2}"],
            # x-j
            "xjmm": ["j-1           "],
            "xjm":  ["j-\\frac{1}{2}"],
            "xjc":  ["j             "],
            "xjp":  ["j+\\frac{1}{2}"],
            "xjpp": ["j+1           "],
            # x-k
            "xkmm": ["k-1           "],
            "xkm":  ["k-\\frac{1}{2}"],
            "xkc":  ["k             "],
            "xkp":  ["k+\\frac{1}{2}"],
            "xkpp": ["k+1           "],
            # y-i
            "yimm": ["i-1           "],
            "yim":  ["i-\\frac{1}{2}"],
            "yic":  ["i             "],
            "yip":  ["i+\\frac{1}{2}"],
            "yipp": ["i+1           "],
            # y-j
            "yjmm": ["j-\\frac{1}{2}"],
            "yjm":  ["j             "],
            "yjc":  ["j+\\frac{1}{2}"],
            "yjp":  ["j+1           "],
            "yjpp": ["j+\\frac{3}{2}"],
            # y-k
            "ykmm": ["k-1           "],
            "ykm":  ["k-\\frac{1}{2}"],
            "ykc":  ["k             "],
            "ykp":  ["k+\\frac{1}{2}"],
            "ykpp": ["k+1           "],
            # z-i
            "zimm": ["i-1           "],
            "zim":  ["i-\\frac{1}{2}"],
            "zic":  ["i             "],
            "zip":  ["i+\\frac{1}{2}"],
            "zipp": ["i+1           "],
            # z-j
            "zjmm": ["j-1           "],
            "zjm":  ["j-\\frac{1}{2}"],
            "zjc":  ["j             "],
            "zjp":  ["j+\\frac{1}{2}"],
            "zjpp": ["j+1           "],
            # z-k
            "zkmm": ["k-\\frac{1}{2}"],
            "zkm":  ["k             "],
            "zkc":  ["k+\\frac{1}{2}"],
            "zkp":  ["k+1           "],
            "zkpp": ["k+\\frac{3}{2}"],
            #
            "dintrpa": ["{\\overline {#1}^{\\left( A,#2 \\right)}}", 2], # discrete arithmetic average
            "dintrpv": ["{\\overline {#1}^{\\left( V,#2 \\right)}}", 2], # discrete volume average
            "dintrpu": ["{\\overline {#1}^{\\left( U,#2 \\right)}}", 2], # discrete average, unknown
        }
    }
}

