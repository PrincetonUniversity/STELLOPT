# GitHub Pages
This is the branch to automatically generate GitHub Pages at https://princetonuniversity.github.io/STELLOPT/.

1. The basic syntax is markdown [(a quick cheatsheet for markdown)](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet).

2. [MathJax]() is configured. Put `{% include head.html %}` at the head of each file and you can type equations by as following
```
{% include head.html %}

inline
$$ J \times B = \nabla B $$
or displayed
\[ J \times B = \nabla B \]
\$$ J \times B = \nabla B $$
```
And you should get the equations printed out on the [GitHub Pages](https://princetonuniversity.github.io/STELLOPT/README) (might not work in the preview).

inline
$$ J \times B = \nabla B $$
or displayed
\[ J \times B = \nabla B \]
\$$ J \times B = \nabla B $$

3. Table of contents can be achieved by inserting links. A easier tool is to use [gh-md-toc](https://github.com/ekalinin/github-markdown-toc).
```markdown
BEAMS3D
=======

* [Theory](#theory)

* [Compilation](#compilation)

* [Input Data Format](#input-data-format)

* [Execution](#execution)

* [Output Data Format](#output-data-format)

* [Visualization](#visualization)

* [Tutorials](#tutorials)
```

4. [script.sh](script.sh) is a batchscript for converting **Creole** format to **markdown** using pandoc. You can modify it doing converts for other formats.
