window.MathJax = {
  loader: {load: ['tex-mml-chtml','[tex]/physics']},
  tex: {
    inlineMath: [["\\(", "\\)"]],
    displayMath: [['$$', '$$'],['$$\n', '\n$$'], ['$$\r\n', '\r\n$$'], ['\\[', '\\]']],
    packages: {'[+]': ['ams','physics']},
    tags: 'ams',
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*",
    processHtmlClass: "arithmatex"
  }
};

document$.subscribe(() => { 
  MathJax.startup.output.clearCache()
  MathJax.typesetClear()
  MathJax.texReset()
  MathJax.typesetPromise()
})