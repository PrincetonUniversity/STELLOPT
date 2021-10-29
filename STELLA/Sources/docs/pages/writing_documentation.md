---
title: Writing Documentation
---

[TOC]

stella uses [FORD][ford] to automatically build this online documentation. There are two
sorts of documentation that get built: in-source and out-of-source. This page describes
how to write both sorts. The website will get built automatically when your changes are
merged into `master`.

# In-source documentation

Special Fortran comments of the form `!>` are used to document procedures, variables and
programs using [FORD][ford]. These comments are sometimes called "docstrings". FORD can
understand these comments before, after or in-line with the thing to be documented, but
the preferred style in stella is to put the comment before the entity, like in this example:

``` fortran
!> This is module-level documentation, describing the overall purpose of this module.
!> Documentation comments can be over several lines.
module my_module
    implicit none

    !> Docstrings can go on derived types
    type :: my_type
        !> Docstrings can include in-line LaTeX like this: \(R\)
        real :: R_major
        !> Or as a displayed equation:
        !>   $$\frac{\partial R}{\partial \psi}$$
        real :: R_major_prime
        !> Code-formatting uses `backticks`
        integer :: n_R
    end type my_type

contains

    !> We can document functions/subroutines
    real function gradient(R, psi)
        !> In order to separate the docstring for each argument...
        real, intent(in) :: R
        !> ...it's best to have each argument on its own line
        real, intent(in) :: psi
    end function gradient

end module my_module
```

The [FORD wiki][ford_wiki] has more documentation on how to write these docstrings. Please
note that while FORD accepts `!<` comments _after_ or in-line with the thing to be
documented, the stella style is to stick to `!>` _before_ and on a separate line.

# Out-of-source documentation

As well as the code documentation, we also have some extra documentation. The source for
these pages in still kept in the stella git repository, and are built into the website at the
same time. These extra pages are written in [Markdown][markdown-guide], with some
extensions (see [the Python markdown implementation][python-markdown]). There's a short
[cheat sheet](#cheat-sheet) at the bottom of this page.

## Add a new page

All the out-of-source documentation is under `docs/pages`. FORD converts the directory
hierarchy in to a hierarchy of HTML pages. To add a new page, simply create a new file
under `docs/pages`, and use the file extension `.md`. **Your file must contain the
following at the very top**:

``` markdown
---
title: Page title
---
```

Without this metadata section, FORD will not parse the file as part of the documentation.

In this metadata section, you can also have `author` and `date` items.

Another useful feature, is put `[TOC]` on its own line after the metadata section. This
produces a hyperlinked table of contents.

See the [FORD wiki][ford_pages] for a more detailed description.

## Converting LaTeX to Markdown

The easiest way to convert LaTeX to Markdown is to use [Pandoc][pandoc]. This nifty tool
understands tons of text formats and can convert between them easily, and will get the
vast majority of the heavy lifting done for you. Due to the complexities of LaTeX and the
simplicity of Markdown, it may require some manual tidying up after the initial
conversion.

To get started, run Pandoc:

``` shell
$ pandoc --standalone --from=latex --to=gfm my_docs.tex --output docs/pages/my_docs.md
```

FORD is a bit fussier about the metadata section at the top, compared to what Pandoc
produces, so you may need to manually adjust it.

Acceptable:

``` markdown
---
title: On the Electrodynamics of Moving Bodies
author: A. Einstein
---
```

Unacceptable:

``` markdown
---
title: Does the Inertia of a Body Depend upon its Energy Content?
author: 
- 'A. Einstein'
---
```

which is a possible output from Pandoc.

FORD can render LaTeX included in the markdown, with just a couple of gotchas. The most
important one is that in-line maths **must** use `\(...\)` rather than `$...$`. Displayed
equations can be written between `$$...$$`, but note that this does not number the
resulting equation.

Normal `\begin{equation} ... \end{equation}` environments can be used to get numbered
equations, along with `\label{eq:something}` and `\eqref{eq:something}` to refer to them.

This example:

```
To obtain the distribution function at the next time step, \(g^{n+1}\), we could combine
these equations

$$A g^{n+1} + B g^{n} = DF^{-1}Gg^{n+1} + E\phi^{n},$$
```

is rendered as:

> To obtain the distribution function at the next time step, \(g^{n+1}\), we could combine
> these equations
>
> $$A g^{n+1} + B g^{n} = DF^{-1}Gg^{n+1} + E\phi^{n},$$

while

```
\begin{equation}
F \phi^{n+1} = G g^{n+1}
\label{eq:QN}
\end{equation}

\eqref{eq:QN} is the quasi-neutrality equation
```

is rendered as

> \begin{equation}
> F \phi^{n+1} = G g^{n+1}
> \label{eq:QN}
> \end{equation}
> 
> \eqref{eq:QN} is the quasi-neutrality equation


## Linking to sections

Section titles within a page get converted into HTML "anchors" which can be linked to. The
section names are first converted to lowercase and spaces replaced with hyphens. To link
to a section, use the usual link syntax and add `#` in front of the converted section
name.

This:

``` markdown
Link to [this section](#linking-to-sections)
```

becomes:

> Link to [this section](#linking-to-sections)

## Linking to source documentation

Linking directly to the code documentation is possible using FORD's syntax, which is
described
[here](https://github.com/Fortran-FOSS-Programmers/ford/wiki/Writing-Documentation#links):

This:

``` markdown
The two linear steps \(L\) are performed by the function [[dist_fn:invert_rhs]] in [[dist_fn.fpp(file)]]
```

is rendered as:

> The two linear steps \(L\) are performed by the function [[dist_fn:invert_rhs]] in [[dist_fn.fpp(file)]]

## Building the documentation locally

[FORD][ford] can be easily installed with `pip`:

``` shell
$ pip3 install --user ford
```

After installing FORD, simply run `make doc`. This will build the documentation under
`docs/html`:

``` shell
$ make doc
```

Then open `docs/html/index.html` in your favourite browser.

## Cheat sheet

[See here][markdown-guide] for a comprehensive guide to Markdown.

Here's a quick little cheat sheet:

| Syntax          | Description                |
|-----------------|----------------------------|
| Heading         | `# Top-level`              |
|                 | `## Section`               |
|                 | `### Sub-section`          |
| Bold            | `**bold text**`            |
| Italic          | `*italicised text*`        |
| Code            | `` `code` ``               |
| Blockquote      | `> Block quote`            |
| Ordered lists   | `1. First item`            |
|                 | `2. Second item`           |
|                 | `3. Third item`            |
| Unordered lists | `- First item`             |
|                 | `- Second item`            |
|                 | `- Third item`             |
| Link            | `[title](www.example.com)` |
| Image           | `![alt text](image.jpg)`   |


[ford]: https://github.com/Fortran-FOSS-Programmers/ford
[ford_wiki]: https://github.com/Fortran-FOSS-Programmers/ford/wiki
[ford_pages]: https://github.com/Fortran-FOSS-Programmers/ford/wiki/Writing-Pages
[markdown-guide]: https://www.markdownguide.org
[pandoc]: https://pandoc.org/
[python-markdown]: https://python-markdown.github.io/

<!-- Local Variables: -->
<!-- mode: gfm -->
<!-- fill-column: 90 -->
<!-- End: -->
