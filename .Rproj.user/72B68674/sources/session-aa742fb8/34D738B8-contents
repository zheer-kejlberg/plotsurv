# {publistR}
Zheer Kejlberg Al-Mashhadi

## Welcome to {publistR}

An R package to easily turn your bibliography or a list of DOIs into a
structured and formatted list of publications.

For an example, see
<a href="https://zheer.dk/cv/publications.html">zheer.dk/cv/publications</a>.
Note, .pdf and .docx versions are also available in the sidebar.

## Installation

First, you need to have the `{devtools}` package installed:

``` r
install.packages("devtools")
```

Then, to install this package from GitHub:

``` r
devtools::install_github("zheer-kejlberg/publistR")
library(publistR)
```

And you’re good to go.

## Usage

The documentation should help you get started quickly:

``` r
?publistR
```

The package consists of one function, also named `publistR()`:

``` r
publistR(
  ref_sections,
  title = "",
  subtitle = "",
  author = "",
  prepend = NULL,
  append = NULL,
  author_names = list(list(family = "Family", given = "Given")),
  merge_sections = FALSE,
  merged_title = NULL,
  custom_fonts = NULL,
  custom_csl = NULL,
  bib_file = NULL,
  title_bold = FALSE,
  title_italic = FALSE,
  title_underline = FALSE,
  title_small_caps = FALSE,
  output_format = "docx",
  output_path = getwd(),
  output_filename = "publication_list",
  add_yaml = NULL,
  keep_files = FALSE
)
```

The argument `ref_sections` is the only required argument.

### Sections

The `ref_sections` argument is what creates the contents of the final
document. Its input must be a list of named lists. Each of these named
lists corresponds to one separate section in the final document. Each
named list must have two elements; a “title” element with one character
value, and a “DOIs” element with a (vector of) character value(s):

``` r
ref_sections = list(
  list(title = "First authorships",
       DOIs = c("https://doi.org/10.3389/fendo.2022.882998",
                "https://doi.org/10.3389/fendo.2022.861422"
       )),
  list(title = "Co-authorships",
       DOIs = c("https://doi.org/10.1016/j.jacc.2020.11.059",
                "https://doi.org/10.1016/j.jacc.2019.06.057"
       ))
)
```

Note, the `DOIs` inputs can be either actual DOI addresses (in short
form or URL form) or references to a local .bib file (more info on this
below).

### Author names in bold font

The `author_names` argument takes a list of named lists as input. Each
named list represents an author name to be highlighted in all citations
in the document. E.g.;

``` r
author_names = list(
  list(family = "Al-Mashhadi", given = "Zheer"),
  list(family = "Al-Mashhadi", given = "Zheer Kejlberg")
)
```

### Using .bib file instead of (or in addition to) DOIs

The `bib_file` argument allows you to specify the filepath to a .bib
file. Any inputs to a `DOIs` key which does not have the required
formatting of a DOI is assumed instead to be a citation key within the
specified .bib file.

E.g.:

``` r
publistR(
  author_names = list(
    list(family = "Al-Mashhadi", given = "Z. K.")
  ),
  ref_sections = list(
    list(title = "Manuscripts", 
         DOIs = c("https://doi.org/10.3389/fendo.2022.882998",
                  "al-mashhadi2022"
         ))
  ),
  bib_file = "references.bib"
)
```

### Other

- `title =`: This takes a single character value for the page title.
  Defaults to blank.
- `subtitle =`: This takes a single character value for the page
  subtitle. Defaults to blank.
- `author =`: This takes a single character value for the page author
  Defaults to blank.
- `prepend =`: Takes a vector of character values, each element
  representing one line in Quarto format, to go before all the
  bibliography sections. Defaults to blank. E.g.:

``` r
prepend = c(
  "\\",
  "# Extra section",
  "\\",
  "Some text",
  "\\"
)
```

- `append =`: Same as prepend but is inserted at the bottom of the
  document. Defaults to blank.
- `merge_sections =`: Convenience argument which allows you to merge all
  sections specified in the `ref_sections` argument to one combined
  section – just to save you some time. A `merged_title` will be
  required.
- `merged_title =`: A title for the merged section from
  `merge_sections`. Ignored if sections aren’t merged.
- `custom_fonts =`: Filepath to a .docx file with custom fonts named
  “font-ref-doc.docx”. Will be used to set the fonts of any .docx output
  files. Defaults to publistR’s internal template.
- `custom_csl =`: Filepath to a .csl (citation style language) file for
  the formatting of the citations. Note, all formatting must be done
  within the .csl, and inputs to the formatting arguments below
  (`title_bold`, `title_italic` etc.) will be ignored. Defaults to an
  apa-based csl.
- `title_bold =`: Boolen (TRUE/FALSE). Make the titles of citations
  bold. Defaults to FALSE.
- `title_italic =`: Boolen (TRUE/FALSE). Make the titles of citations
  italic. Defaults to FALSE.
- `title_underline =`: Boolen (TRUE/FALSE). Make the titles of citations
  underlined. Defaults to FALSE.
- `title_small_caps =`: Boolen (TRUE/FALSE). Make the titles of
  citations small-caps formatted. Defaults to FALSE.
- `output_format =`: Choose the filetype for the document to be
  generated. Options are “pdf”, “html”, “docx” or “all”. Defaults to
  “docx”.
- `output_path =`: Filepath, where to save the generated document.
  Defaults to current working directory.
- `output_filename =`: Choose the name of the generated document
  (*without* filename extension, e.g., “docx”). Defaults to
  “publication_list”.
- `add_yaml =`: Advanced option. Supply a list of lists in YAML format
  to be added to the intermediately generated Quarto document. Top-level
  keys *cannot* be one of “title”, “subtitle”, “author”, “format”,
  “editor”, “section-bibs-bibliography”, “csl”, “citeproc”, or
  “filters”. To alter these, see the final argument `keep_files`.
- `keep_files =`: Boolean (TRUE/FALSE). If set to TRUE, the
  intermediately generated files are not deleted. Can be used for
  debugging purposes or if more fine-grained control is needed (e.g., of
  the YAML section of the .qmd document). Defaults to FALSE.

Enjoy.

## Acknowledgements

This package draws on open source work by the following:

Brenton Wiernik’s APA .csl:
https://github.com/citation-style-language/styles/blob/master/apa-cv.csl

Shafee’s bold-author.lua file from the following SO thread:
https://stackoverflow.com/questions/76394078/format-specific-authors-with-bold-font-in-bibliography-created-with-quarto/76429867#76429867

Matt Capaldi’s hide-me.lua file

Jesse Rosenthal and Albert Krewinkel’s section-bibliographies.lua
