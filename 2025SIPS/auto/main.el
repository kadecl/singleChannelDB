(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("graphicx" "dvipdfmx") ("xcolor" "dvipdfmx")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "graphicx"
    "xcolor"
    "amsmath"
    "amsthm"
    "amssymb"
    "latexsym"
    "comment"
    "lipsum"
    "ulem"
    "cancel"
    "bm"
    "mathtools"
    "multirow"
    "cuted"
    "algorithm"
    "algorithmic"
    "booktabs"
    "lscape")
   (TeX-add-symbols
    "argmax"
    "argmin"
    "minimize")
   (LaTeX-add-environments
    "thm"
    "lem"
    "prop"
    "cor"
    "dfn"
    "rem"
    "prob"))
 :latex)

