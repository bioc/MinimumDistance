(TeX-add-style-hook "MinimumDistance"
 (lambda ()
    (TeX-add-symbols
     '("Rclass" 1)
     '("Rpackage" 1)
     '("Robject" 1)
     '("Rcode" 1)
     '("Rmethod" 1)
     '("Rfunction" 1)
     "R"
     "md"
     "logRratio"
     "baf"
     "bafs"
     "tsl"
     "ts")
    (TeX-run-style-hooks
     "geometry"
     "url"
     "natbib"
     "graphicx"
     "latex2e"
     "art10"
     "article")))

