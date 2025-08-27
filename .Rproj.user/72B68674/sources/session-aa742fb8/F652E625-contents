#' publistR
#'
#' @name publistR
#' @rdname publistR
#' @author Zheer Kejlberg Al-Mashhadi
#' @description Create publication list
#' @export
#' @usage publistR(
#'          ref_sections,
#'          title = "",
#'          subtitle = "",
#'          author = "",
#'          prepend = NULL,
#'          append = NULL,
#'          author_names = list(list(family = "Family", given = "Given")),
#'          merge_sections = FALSE,
#'          merged_title = NULL,
#'          custom_fonts = NULL,
#'          custom_csl = NULL,
#'          bib_file = NULL,
#'          title_bold = FALSE,
#'          title_italic = FALSE,
#'          title_underline = FALSE,
#'          title_small_caps = FALSE,
#'          output_format = "docx",
#'          output_path = getwd(),
#'          output_filename = "publication_list",
#'          add_yaml = NULL,
#'          keep_files = FALSE
#'        )
#' @return NULL

#' @param ref_sections a list of lists containing a section title and DOIs for each section of the reference paper. Each embedded list must have two keys, "title =" and "DOIs =".
#' @param title (optional) a character value for the document title
#' @param subtitle (optional) a character value for the document subtitle
#' @param author (optional) a character value for the document author
#' @param prepend (optional) a (vector of) character string(s) to go in the bottom of the final output document, each element in the vector representing one (Quarto formatted) line.
#' @param append (optional) a (vector of) character string(s) to go in the top of the final output document, each element in the vector representing one (Quarto formatted) line.
#' @param author_names (optional) a list of lists containing author names to be highlighted. Each embedded list must have to keys "family =" and "given =" for sur- and firstname, respectively.
#' @param merge_sections (optional) a Boolean (defaults to FALSE) determining whether to merge all the supplied sections into one section with a common section title. For ease of use.
#' @param merged_title (optional) a section title for the merged section IF merge_sections is set to TRUE.
#' @param custom_fonts (optional) File path to a custom Quarto Word Template (for adjusting the output fonts). The file must be named "font-ref-doc.docx". The fonts to set are Heading 1, Hyperlink and Body Text/First Paragraph. Defaults to 'NULL', using PublistR's internal Word Template.
#' @param custom_csl (optional) File path to a custom .csl file to change the citation style. Note: author name bolding is only guaranteed to work with publistR's internal .csl file.
#' @param bib_file (optional) File path to a provided .bib file from which one wishes to reference certain elements. Can be used in conjunction with DOIs.
#' @param title_bold (optional) Boolean, defaults to FALSE; bold reference title. Note: is not used if a custom_csl is set.
#' @param title_italic (optional) Boolean, defaults to FALSE; italic reference title. Note: is not used if a custom_csl is set.
#' @param title_underline (optional) Boolean, defaults to FALSE; underlines reference title. Note: is not used if a custom_csl is set.
#' @param title_small_caps (optional) Boolean, defaults to FALSE; reference title in small caps. Note: is not used if a custom_csl is set.
#' @param output_format (optional) a string, can take either "pdf", "docx", "html" or "all"; defaults to "docx"
#' @param output_path (optional) the path where the output document should be saved, defaults to the working directory
#' @param output_filename (optional) the filename for the final document; defaults to "publication_list"
#' @param add_yaml (optional, advanced) a nested list that can be turned into a YAML object (see yaml::as.yaml()) with additional metadata parameters. Some keys are already taken and can not be customised.
#' @param keep_files (optional) Boolean, determines whether temporary files created in the process should be deleted
#' @examples
#'   \dontrun{
#'   publistR(
#'     ref_sections = list(
#'       list(title = "First authorships",
#'            DOIs = c("https://doi.org/10.3389/fendo.2022.882998",
#'                     "https://doi.org/10.3389/fendo.2022.861422"
#'            )),
#'       list(title = "Co-authorships",
#'            DOIs = c("https://doi.org/10.1016/j.jacc.2020.11.059",
#'                     "https://doi.org/10.1016/j.jacc.2019.06.057"
#'            ))
#'     ),
#'     author_names = list(
#'       list(family = "Al-Mashhadi", given = "Z."),
#'       list(family = "Al-Mashhadi", given = "Z. K.")
#'     ),
#'     title_italic = T,
#'     output_format = "pdf"
#'   )
#'   }

#### publistR():

publistR <- function(ref_sections,
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
) {
  #### ERROR HANDLING ####
  if (output_format != "pdf" & output_format != "docx" & output_format != "html" & output_format != "all") {
    stop("Please select an approved output format")
  }
  if (is.null(ref_sections)) {
    stop("No DOIs supplied")
  }
  if (merge_sections == TRUE & is.null(merged_title)) {
    stop("merge_sections argument set to TRUE but no input supplied to the merged_title argument. The merged section must have a title.")
  }
  # check keys for author_names
  check_author_name_keys <- function(input) {
    input == c("family", "given")
  }
  if (!all(unlist(lapply(lapply(author_names, names), check_author_name_keys)))) {
    stop("Author names keys not specified correctly. Make sure there are only key pairs of 'family' and 'given'.")
  }
  # add comma to given names
  for (name_i in 1:length(author_names)) {
    author_names[[name_i]][["given"]] <- paste0(author_names[[name_i]][["given"]],",")
  }
  # check keys for ref_sections
  check_ref_sections_keys <- function(input) {
    input == c("title", "DOIs")
  }
  if (!all(unlist(lapply(lapply(ref_sections, names), check_ref_sections_keys)))) {
    stop("Ref_section keys not specified correctly. Make sure there are only key pairs of 'title' and 'DOIs'.")
  }
  # Check if custom_fonts word template exists
  if (!is.null(custom_fonts)) {
    if (!file.exists(custom_fonts)) {
      stop(paste0("Check that the required Word Template exists at ", custom_fonts))
    }
  }
  # Check if custom_csl file exists
  if (!is.null(custom_csl)) {
    if (!file.exists(custom_csl)) {
      stop(paste0("Check that the required .csl file exists at ", custom_csl))
    }
  }
  # Check YAML keys and values
  if (length(title) > 1) { stop("Please specify a single character value in `title =`") }
  if (length(subtitle) > 1) { stop("Please specify a single character value in `subtitle =`") }
  if (length(author) > 1) { stop("Please specify a single character value in `author =`") }
  for (yaml_line in add_yaml) {
    if (names(yaml_line) %in% c("title", "subtitle", "author", "format", "editor", "section-bibs-bibliography", "csl", "citeproc", "filters")) {
      stop(paste0("'", names(yaml_line), "' is not an approved custom YAML key."))
    }
  }


  #### CREATE NECESSARY FILES ####
  # get filepaths
  if (!is.null(custom_csl)) {
    csl_path <- custom_csl
  } else {
    csl_path <- system.file("publistR.csl", package = "publistR")
  }

  bold_author_path <- system.file("bold-author.lua", package = "publistR")
  hide_me_path <- system.file("hide-me.lua", package = "publistR")
  if (!is.null(custom_fonts)) {
    font_doc_path <- custom_fonts
  } else {
    font_doc_path <- system.file("font-ref-doc.docx", package = "publistR")
  }
  section_bibliographies_path <- system.file("section-bibliographies.lua", package = "publistR")

  # location to copy to
  to_path <- paste0(getwd(), "/publistR_temp")

  # create temp directory
  dir.create(to_path)

  # copy files
  file.copy(bold_author_path, to_path)
  file.copy(hide_me_path, to_path)
  file.copy(font_doc_path, to_path)
  file.copy(section_bibliographies_path, to_path)

  #### MODIFY CSL ####
  input_formatting <- function(input) {
    for (line in c(534,540,546,549,554)) {
      csl[line] <- gsub("/>", paste0(" ", input," />"), csl[line])
    }
    return(csl)
  }
  csl <- readLines(csl_path)
  if (is.null(custom_csl)) {
    if (title_bold) { csl <- input_formatting("font-weight=\"bold\"") }
    if (title_italic) { csl <- input_formatting("font-style=\"italic\"") }
    if (title_underline) { csl <- input_formatting("text-decoration=\"underline\"") }
    if (title_small_caps) { csl <- input_formatting("font-variant=\"small-caps\"") }
  }
  writeLines(csl, paste0(to_path, "/publistR.csl"))


  #### MODIFY YAML HEADER ####
  # Read YAML file
  yaml_path <- system.file("yaml.yaml", package = "publistR")
  yaml_data <- yaml::yaml.load_file(yaml_path)

  # add author names
  if (!is.null(author_names)) {
    yaml_data[["bold-auth-name"]] <- author_names
  }
  yaml_data[["title"]] <- title
  yaml_data[["subtitle"]] <- subtitle
  yaml_data[["author"]] <- author

  # Add custom YAML
  if (!is.null(add_yaml)) {
    for (yaml_line in add_yaml) {
      yaml_data[[names(yaml_line)]] <- unname(yaml_line)
    }
  }

  # fix the citeproc value
  yaml_data <- yaml::as.yaml(yaml_data, handlers = list(logical=yaml::verbatim_logical))

  # Write the modified YAML to a new file in the temp directory
  temp_path <- paste0(getwd(), "/publistR_temp/yaml.txt")
  writeLines(yaml_data, temp_path)

  #### ADD REFERENCES ####
  # Load the supplied .bib-file if any
  if (!is.null(bib_file)) {
    custom_bib <- paste0(readLines(bib_file), collapse="")
  }
  #Get paper info by DOI
  message("Retrieving publication info from supplied DOIs.")
  extract_dois <- function(ref_section) {
    output <- list()
    for (i in 1:length(ref_section$DOIs)) {
      DOI <- ref_section$DOIs[i]
      if(grepl("10.[0-9]+/.+", DOI) | is.null(bib_file)) {
        output[i] <- rcrossref::cr_cn(DOI)
      } else {
        # get the ref
        pattern <- paste0("@[^@]*?\\{",DOI,",.*?(?=@)")
        output[i] <- unlist(stringr::str_extract_all(custom_bib, pattern))
      }
    }
    return(output)

  }
  bibtex <- lapply(ref_sections,extract_dois)
  message("DOIs retrieved.")
  # Rename shorthand/keys

  rename_doi_key <- function(dois) {
    for (i in 1:length(dois)) {
      doi <- stringr::str_extract(dois[i], "DOI=\\{.+?\\}")
      doi <- gsub("DOI=\\{","",doi)
      doi <- gsub("\\}","",doi)
      dois[i] <- gsub("\\{.+, title=\\{", paste0("\\{",doi,", title=\\{"), dois[i])
    }
    dois <- as.list(dois)
  }
  bibtex <- lapply(bibtex,rename_doi_key)

  # Create sections
  get_section_titles <- function(ref_sections) {
    ref_sections$title
  }
  ref_section_titles <- lapply(ref_sections,get_section_titles)

  if (merge_sections == FALSE) {
    sections <- c()
  } else {
    sections <- c("", paste0("# ", merged_title, "  {.sectionbibliography}"), "")
  }

  for (i in 1:length(ref_section_titles)) {

    #get_shorthand <- function(DOI) { sub(", (title|type)[ ]*=.*$", "", sub("^[ ]?@.+?\\{", "", DOI)) }
    get_shorthand <- function(DOI) { sub("^[ ]*?@.+?\\{", "", DOI) }

    if (merge_sections == FALSE) {
      sections <- c(sections, "", paste0("# ", unlist(ref_section_titles[i]), "  {.sectionbibliography}"), "","::: hide-me", paste0("@", get_shorthand(unlist(bibtex[i]))),":::", "", "::: {.sectionrefs}", ":::", "")
    } else {
      sections <- c(sections, "::: hide-me", paste0("@", get_shorthand(unlist(bibtex[i]))),":::", "", "::: {.sectionrefs}", ":::", "")
    }

  }

  # Create .bib file
  bib <- unlist(bibtex)
  bib <- stringr::str_replace_all(bib,"‐","-")
  writeLines(bib, paste0(getwd(), "/publistR_temp/references.bib"))

  #### CREATE .QMD FILE ####
  yaml_section <- readLines(paste0(getwd(), "/publistR_temp/yaml.txt"))
  yaml_section <- c("---", yaml_section, "---","")

  if (!is.null(prepend)) { yaml_section <- c(yaml_section, prepend) }
  qmd_file <- c(yaml_section, "", sections, "")
  if (!is.null(append)) { qmd_file <- c(qmd_file, append) }

  qmd_filename <- paste0(output_filename,".qmd")
  writeLines(qmd_file,
             paste0(getwd(), paste0("/",qmd_filename))
             )

  #### KNIT .QMD FILE ####
  saved_wd <- getwd() # save current wd first
  setwd(output_path) # set new wd for the output
  if (output_format != "all") {
    output_filename <- paste0(output_filename, ".", output_format)
  }
  quarto::quarto_render(input = paste0(saved_wd, paste0("/",qmd_filename)),
                        output_format = output_format)

  #### DELETE HELPER FILES ####
  setwd(saved_wd) # return to saved wd to clean up

  if (!keep_files) {
    file.remove(qmd_filename)
    unlink("publistR_temp", recursive = TRUE)
  }
}
