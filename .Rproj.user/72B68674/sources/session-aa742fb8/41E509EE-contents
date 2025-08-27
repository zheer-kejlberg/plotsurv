#' plotsurv
#'
#' @name plotsurv
#' @rdname plotsurv
#' @author Zheer Kejlberg Al-Mashhadi
#' @description Create nice plots for Kaplan-Meier/Aalen-Johansen models fitted using survival::survfit().
#' @export
#' @usage plotsurv(survfit_obj,
#'          include_surv = TRUE,
#'          include_CIFs = NULL,
#'          display_event = "all",
#'          title = "Cumulative incidence and survival",
#'          subtitle = "",
#'          x_lab = "Time",
#'          y_lab = "Cumulative probability of event",
#'          color_lab = "Group and event type",
#'          fill_lab = "Group and event type",
#'          line_lab = "Group and event type",
#'          group_labels = NULL,
#'          colors = NULL,
#'          linetypes = NULL
#'        )

#' @return ggplot object

#' @param survfit_obj the output of a survift() call
#' @param include_surv (optional) whether to overlay the Kaplan-Meier survival curve
#@param include_CIFs (optional)
#' @param conf.int (optional) whether to add confidence bands (defaults to TRUE)
#' @param ticks (optional) whether to add censoring ticks (defaults to TRUE)
#' @param ticksize (optional) size of censoring ticks (defaults to 3)
#' @param tickalpha (optional) colour alpha-level of censoring ticks (defaults to 0.8)
#' @param display_event (optional) a vector of event values to be displayed (defaults to "all"),
#' @param title (optional) The plot title, defaults to "Cumulative incidence and survival"
#' @param subtitle (optional) The plot subtitle, defaults to empty string ""
#' @param x_lab (optional) X-axis label, defaults to "Time"
#' @param y_lab (optional) Y-axis label, defaults to "Cumulative probability of event"
#' @param color_lab (optional) Label for color-codes, defaults to "Group and event type" (if same as fill_lab and line_lab, one combined legend will apear)
#' @param fill_lab (optional) Label for color-fills, defaults to "Group and event type" (if same as color_lab and line_lab, one combined legend will apear)
#' @param line_lab (optional) Label for lines, defaults to "Group and event type" (if same as color_lab and fill_lab, one combined legend will apear)
#' @param group_labels (optional) Labels for each subgroup (defaults to automated combination of RHS group and LHS event type), must be of same length as no. of groups
#' @param colors (optional) Colors for each subgroup, must be of same length as no. of groups
#' @param linetypes (optional) Linetypes for each subgroup, must be of same length as no. of groups
#'
#'
#' @examples
#'   \dontrun{
#'   plotsurv(
#'     data = fit_obj,
#'     include_surv = FALSE,
#'     title = "Cumulative incidence of metastatic breast cancer recurrence",
#'     subtitle = "Between 2015 and 2023",
#'     colors = c("#3A488A","#BD5630")
#'   )
#'   }

#### plotsurv():

plotsurv <- function(survfit_obj, # the output of a call to survival::survfit()
                     include_surv = TRUE, # overlay the Kaplan-Meier survival curve
                     #include_CIFs = NULL, # include the cause-specific cumulative incidence curves (with one outcome type, this is 1 - S(t); otherwise, it is estimated by survfit() using the Aalen-Johansen estimator)
                     conf.int = TRUE, #
                     ticks = TRUE,
                     ticksize = 3,
                     tickalpha = 0.8,
                     display_event = "all", #
                     title = "Cumulative incidence and survival", # Plot title
                     subtitle = "", # Plot subtitle
                     x_lab = "Time", # X-axis label
                     y_lab = "Cumulative probability of event", # Y-axis label
                     color_lab = "Group and event type", # Label for colors
                     fill_lab = "Group and event type", # Label for fillings
                     line_lab = "Group and event type", # Label for linetypes
                     group_labels = NULL, # Labels to rename each group
                     colors = NULL, # Colors
                     linetypes = NULL) # Linetypes
{

  #print(display_event)
  #print(survfit_obj$states)
  if (identical(display_event, "all")) {display_event <- survfit_obj$states} else {display_event <- unique(c("(s0)", display_event))}
  if (! all(display_event %in% survfit_obj$states)) {
    stop("Some values of display_event were not in the event arg of your Surv() call (note: the lowest factor, corresponding to the censoring state, must not be included.)")
  }

  censoring_times <- survfit_obj$time[survfit_obj$n.censor[,"(s0)"] == 1]

  refactor_survfits <- function(survfit_obj) {
    event_types <- survfit_obj$states

    create_dfs <- function(event_type, survfit_obj) {
      df <- data.frame(
        time = survfit_obj$time,
        est = survfit_obj$pstate[,event_type],
        upper = survfit_obj$upper[,event_type],
        lower = survfit_obj$lower[,event_type],
        strata = paste0(rep(names(survfit_obj$strata), survfit_obj$strata), ", event = ", event_type)
      )
    }

    df_list <- lapply(event_types, create_dfs, survfit_obj)
    names(df_list) <- event_types

    return(df_list)
  }

  dfs <- refactor_survfits(survfit_obj)
  #return(dfs)
  dfs <- dfs[display_event]

  plot_obj <- ggplot2::ggplot() +
    labs(
      y = y_lab,
      x = x_lab,
      title = title,
      subtitle = subtitle,
      color = color_lab,
      fill = fill_lab,
      linetype = line_lab
    )



  if (is.null(linetypes)) {
    linetypes <- rep("solid", length(names(dfs)) * length(unique(dfs[[1]]$strata)))
  }

  if (!(is.null(group_labels))) {
    if (is.null(colors)) {  colors <- scales::hue_pal()(length(group_labels))  }
    plot_obj <- plot_obj +
      scale_color_manual(
        labels = group_labels,
        values = colors
      ) +
      scale_fill_manual(
        labels = group_labels,
        values = colors
      ) +
      scale_linetype_manual(
        labels = group_labels,
        values = linetypes
      )
  } else {
    if (!is.null(colors)) {
      plot_obj <- plot_obj +
        scale_color_manual(
          values = colors
        ) +
        scale_fill_manual(
          values = colors
        ) +
        scale_linetype_manual(
          values = linetypes
        )
    }
  }

  add_to_plot <- function(df, event_type) {
    if (event_type != "(s0)" | include_surv) {
      plot_obj <- plot_obj +
        geom_line(data = df,
                  aes(x = time, y = est, color = strata, linetype = strata))
      if (conf.int) {
        plot_obj <- plot_obj +
          geom_ribbon(data = df,
                      aes(x = time, ymin = lower, ymax = upper, fill = strata),
                      alpha = 0.5, color = NA)
      }
      if (ticks) {
        plot_obj <- plot_obj +
          geom_point(data = subset(df, time %in% censoring_times),
                     aes(x = time, y = est, color = strata),
                     alpha=tickalpha,
                     shape="|",
                     size = ticksize)
      }

    }
    return(plot_obj)
  }

  for (event_type in names(dfs)) {
    plot_obj <- add_to_plot(dfs[[event_type]], event_type)
  }

  return(plot_obj)
}
