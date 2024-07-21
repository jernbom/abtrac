
#' Compute distance matrix from antibody trajectories
#'
#'`ab_tr_dist()` computes a distance matrix using the `cosine*euclidean` distance metric.
#'
#' @param trajectories A tibble with columns `trajectory`, `t`, and `fc`. See `abtrac::trajectories` for an example.
#'
#' @return A distance matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' ab_tr_dist(trajectories = trajectories)
#' }
ab_tr_dist <- function(trajectories) {
  trs <-
    trajectories %>%
    dplyr::select(.data$trajectory, .data$t, .data$fc) %>%
    tidyr::pivot_wider(names_from = .data$t, values_from = .data$fc) %>%
    tibble::column_to_rownames("trajectory") %>%
    as.matrix()

  d_euc <- proxy::dist(trs)
  d_cos <- proxy::dist(trs, method = "cosine", FUN = function(x) 1-x)
  dist_mat <- d_cos*d_euc

  return(dist_mat)
}

#' Cluster antibody trajectories
#'
#' `ab_tr_clust()` takes a distance matrix of antibody trajectories (typically the output from `ab_tr_dist()`) and assigns trajectories to `k` clusters across given `ks` using the PAM (Partitioning Around Medioids) algorithm from `cluster::pam()`.
#'
#' @param dist_mat A distance matrix.
#' @param ks A numeric vector of `k`s to iterate across.
#' @param verbose A logical specifying if progress should be printed.
#' @param ... Additional arguments to be passed to `cluster::pam()`.
#'
#' @return A tibble of `pam` objects.
#' @export
#'
#' @examples
#' \dontrun{
#' ab_tr_clust(dist_mat = dist_mat, ks = 1:20, verbose = TRUE, variant = "faster", nstart = 10)
#' }
ab_tr_clust <- function(dist_mat, ks, verbose = TRUE, ...) {
  t_start <- Sys.time()
  pam_objects <- list()

  for (k in ks) {
    t_0 <- Sys.time()
    if (verbose) cat(stringr::str_c(t_0, "; k = ", k, "\n"))
    pam_objects[[as.character(k)]] <-
      cluster::pam(dist_mat, k = k, diss = TRUE, ...)
    t_1 <- Sys.time()
    if (verbose) print(t_1-t_0)
  }
  t_end <- Sys.time()
  pam_objects <- tibble::enframe(pam_objects, name = "k", value = "pam_object") %>%
    dplyr::mutate(k = as.numeric(.data$k))
  if (verbose) {cat("Clustering complete\n"); print(t_end-t_start)}
  return(pam_objects)
}

#' Silhouette plot
#'
#' Makes a `silhouette_plot()` from a tibble of `pam` objects, typically the output of `ab_tr_clust()`. The silhouette method can be used to dplyr::select the optimal cluster number.
#'
#' @param pam_objects A tibble of `pam` objects with columns `k` (numeric) and `pam_object` (list of `pam` objects). Typically the output of `ab_tr_clust()`.
#' @param nudge_y Nudge cluster number label.
#' @param text_size Size of cluster number label.
#'
#' @return A `gg` object.
#' @export
#'
#' @examples
#' \dontrun{
#' silhouette_plot(pam_objects = pam_objects)
#' }
silhouette_plot <- function(pam_objects, nudge_y = 0.015, text_size = 8) {
  sil <-
    pam_objects %>%
    dplyr::rowwise() %>%
    dplyr::mutate(avg_sil_width = list(.data$pam_object$silinfo$avg.width)) %>%
    dplyr::select(.data$k, .data$avg_sil_width) %>%
    tidyr::unnest(.data$avg_sil_width, keep_empty = TRUE)

  p <-
    sil %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$k, y = .data$avg_sil_width)) +
    ggplot2::geom_line() +
    ggplot2::geom_text(ggplot2::aes(label = .data$k), nudge_y = nudge_y, size = text_size/ggplot2::.pt) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "k", y = "Average silhouette width")

  return(p)
}

#' Extract clusters for a dplyr::selected `k`
#'
#' `extract_clusters()` extracts the cluster assignment from a tibble of `pam` objects and adds them to the original tibble of trajectories.
#'
#' @param trajectories The original tibble of trajectories used to compute the distance matrix (typically using `ab_tr_dist()`). See `abtrac::trajectories` for an example.
#' @param pam_objects A tibble of `pam` objects with columns `k` (numeric) and `pam_object` (list of `pam` objects). Typically the output of `ab_tr_clust()`
#' @param sel_k The dplyr::selected cluster number `k`, possibly dplyr::selected using the silhouette method (see `silhouette_plot()`).
#'
#' @return `trajectories` with an added column of `cluster` numbers.
#' @export
#'
#' @examples
#' \dontrun{
#' extract_clusters(trajectories = trajectories, pam_objects = pam_objects, sel_k = 7)
#' }
extract_clusters <- function(trajectories, pam_objects, sel_k) {
  pam_object <-
    pam_objects %>%
    dplyr::filter(.data$k == sel_k) %>%
    dplyr::pull(.data$pam_object) %>%
    purrr::pluck(1)

  pam_clusters <-
    pam_object$clustering %>%
    tibble::enframe(name = "trajectory", value = "cluster") %>%
    dplyr::mutate(trajectory = as.numeric(.data$trajectory))

  trajectories_clust <-
    trajectories %>%
    dplyr::full_join(pam_clusters, by = c("trajectory"))

  return(trajectories_clust)
}

#' Lineplots of clusters
#'
#' Makes a `cluster_lineplot()` to visualize the obtained clusters of trajectories.
#'
#' @param trajectories_clust Tibble of trajectories and clusters, typically the output of `extract_clusters()`. See `abtrac::trajectories` for an example of all required columns except `cluster`, which is added by `extract_clusters()`.
#' @param alpha A number (in 0:1) defining the transparency of lines.
#' @param ncol A number defining the number of columns for facetting.
#' @param nrow A number defining the number of rows for facetting.
#'
#' @return A `gg` object.
#' @export
#'
#' @examples
#' \dontrun{
#' cluster_lineplot(trajectories_clust = trajectories_clust, alpha = 0.2, ncol = 4)
#' }
cluster_lineplot <- function(trajectories_clust, alpha = 0.2, ncol = NULL, nrow = NULL) {

  cluster_counts <-
    trajectories_clust %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::summarise(n_trajectories = length(unique(.data$trajectory))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(label = stringr::str_c("n(traj)", " = ", c(.data$n_trajectories), collapse = "\n"))

  p <-
    trajectories_clust %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$t, y = .data$log2_fc)) +
    ggplot2::geom_line(ggplot2::aes(group = .data$trajectory), alpha = alpha) +
    ggplot2::facet_wrap(facets = ggplot2::vars(.data$cluster), ncol = ncol, nrow = nrow) +
    ggplot2::scale_y_continuous(breaks = ~ c(seq(0, floor(.x[1]), by = -2), seq(0, ceiling(.x[2]), by = 2)), minor_breaks = NULL) +
    ggplot2::scale_x_continuous(breaks = ~ floor(.x[1]):ceiling(.x[2]), minor_breaks = NULL) +
    ggplot2::labs(x = "Timepoint", y = "log2 foldchange") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank()) +

    ggplot2::geom_text(data = cluster_counts, mapping = ggplot2::aes(label = .data$label, x = min(trajectories_clust$t), y = min(trajectories_clust$log2_fc)), vjust = "inward", hjust = "inward", nudge_x = 0.1, nudge_y = 0.5)

  return(p)
}
