as_dag_coordinates = function(...) {
    coord_dag <- data.frame(rbind(...))
    coord_dag$name = rownames(coord_dag)
    coord_dag
}

dagify = function (..., labels, adjusted=c()) {
    if (missing(labels)) {
        dag = ggdag::dagify(...)
        edges = dagitty::edges(dag)
        vertices = c(edges$v, edges$w)
        labels = setNames(vertices, vertices)
    } else {
        dag = ggdag::dagify(..., labels=labels)
    }
    dagitty::adjustedNodes(dag) = adjusted
    attr(dag, 'labels') = labels
    dag
}

ancestors_except = function(dag, node, except) {
    to_visit = node
    ancestors = c()
    already_visited = c()
    while (length(to_visit) > 0) {
        n = to_visit[1]
        already_visited = c(already_visited, n)
        to_visit = to_visit[-1]
        p = dagitty::parents(dag, n)
        p = p[!(p %in% except)]
        ancestors = c(ancestors, p)
        p = p[!(p %in% c(to_visit, already_visited))]
        to_visit = c(to_visit, p)
    }
    ancestors
}

node_dagitty = function(dag, horizontal_edges_nudge=0, effect='total') {
    nodes = as.data.frame(ggdag::node_status(dag))
    if (!('label' %in% colnames(nodes))) {
        nodes$label = nodes$name
    }

    exposure_ancestors = ancestors_except(dag, dagitty::exposures(dag), except=dagitty::outcomes(dag))
    outcome_ancestors = ancestors_except(dag, dagitty::outcomes(dag), except=dagitty::exposures(dag))

    is_exposure = nodes$name %in% dagitty::exposures(dag)
    is_outcome = nodes$name %in% dagitty::outcomes(dag)
    is_latent = nodes$name %in% dagitty::latents(dag)
    is_adjusted = nodes$name %in% dagitty::adjustedNodes(dag)
    is_exposure_ancestor = nodes$name %in% exposure_ancestors
    is_outcome_ancestor = nodes$name %in% outcome_ancestors

    is_ancestor_of_both = is_exposure_ancestor & is_outcome_ancestor

    nodes$kind = dplyr::case_when(
        is_adjusted == T ~ 'adjusted variable',
        is_exposure == T ~ 'exposure',
        is_outcome == T ~ 'outcome',
        is_latent == T ~ 'unobserved (latent)',
        is_ancestor_of_both == T ~ 'ancestor of exposure and outcome',
        is_exposure_ancestor == T ~ 'ancestor of exposure',
        is_outcome_ancestor == T ~ 'ancestor of outcome',
        .default = 'other variable'
    )

    y_nudge = ifelse(nodes$y == nodes$yend, horizontal_edges_nudge, 0)
    nodes$y_nudged = nodes$y + y_nudge
    nodes$yend_nudged = nodes$yend + y_nudge

    paths = dagitty::paths(dag, Z=dagitty::adjustedNodes(dag))

    nodes$path_kind = apply(nodes, FUN=function(edge) {
        from = edge[['name']]
        to = edge[['to']]
        direction = edge[['direction']]
        if (is.na(to)) {
            return(NA)
        }

        bidirectional_paths = gsub('<|>', '', paths$paths)

        path_fragment = paste(from, to, sep=' - ')
        is_matching_path = grepl(path_fragment, bidirectional_paths)
        if (any(is_matching_path)) {
            matching_path_id = which(is_matching_path)
            is_open = paths$open[matching_path_id]
            if (length(is_open) > 1) {
                matching_path_id = matching_path_id[[ifelse(any(is_open), which(is_open)[[1]], 1)]]
                is_open = any(is_open)
            }
            path = paths$paths[[matching_path_id]]
        } else {
            path_fragment = paste(to, from, sep=' - ')
            is_matching_path = grepl(path_fragment, bidirectional_paths)
            if (any(is_matching_path)) {
                matching_path_id = which(is_matching_path)
                is_open = paths$open[matching_path_id]
                if (length(is_open) > 1) {
                    #stopifnot(all(is_open) || all(!is_open))
                    matching_path_id = matching_path_id[[ifelse(any(is_open), which(is_open)[[1]], 1)]]
                    is_open = any(is_open)
                }
                path = paths$paths[[matching_path_id]]
            } else {
                is_open = dagitty::dconnected(dag, from, to)
                path = paste(from, direction, to)
            }
        }
        fragments = strsplit(path, ' ')[[1]]

        if (is_open) {
            involves_outcome_and_exposure = (
                dagitty::exposures(dag) %in% fragments
                && dagitty::outcomes(dag) %in% fragments
            )
            if (!involves_outcome_and_exposure) {
                # TODOL also fails on Sebastiani
                already_adjusted = dagitty::adjustedNodes(dag)
                required_adjustments = dagitty::adjustmentSets(
                    dag,
                    exposure=dagitty::exposures(dag),
                    outcome=dagitty::outcomes(dag),
                    effect=effect
                )

                is_possibly_biasing = (
                    length(required_adjustments) > 0 &&
                    any(sapply(required_adjustments, function(set) {
                        still_require_adjusting = setdiff(set, already_adjusted)
                        from %in% still_require_adjusting
                    }))
                )

                affects_exposure = dagitty::exposures(dag) %in% fragments && (
                    igraph::distances(
                        igraph::graph_from_data_frame(
                            dagitty::edges(paste('dag{', path, '}'))
                        ),
                        v=to,
                        to=dagitty::exposures(dag),
                        mode='out'
                    ) < Inf
                )

                affects_outcome = dagitty::outcomes(dag) %in% fragments && (
                    igraph::distances(
                        igraph::graph_from_data_frame(
                            dagitty::edges(paste('dag{', path, '}'))
                        ),
                        v=to,
                        to=dagitty::outcomes(dag),
                        mode='out'
                    ) < Inf
                )

                if ((is_possibly_biasing && affects_exposure) || (is_possibly_biasing && affects_outcome)) {
                    return('biasing path')
                }

                # TODO: this fails on Thoemmes example (does not highlight e0)
                return(NA)
            }
            is_causal_path = (
                igraph::distances(
                    igraph::graph_from_data_frame(
                        dagitty::edges(paste('dag{', path, '}'))
                    ),
                    v=dagitty::exposures(dag),
                    to=dagitty::outcomes(dag),
                    mode='out'
                ) < Inf
            )
            required_adjustments = dagitty::adjustmentSets(
                paste('dag{', path, '}'),
                exposure=dagitty::exposures(dag),
                outcome=dagitty::outcomes(dag),
                effect=effect
            )
            already_adjusted = dagitty::adjustedNodes(dag)
            is_well_adjusted = any(sapply(required_adjustments, function(set) {
                still_require_adjusting = setdiff(set, already_adjusted)
                length(still_require_adjusting) == 0
            }))
            is_biasing_path = !is_well_adjusted
            # when a path is both causal and biasing, show it as biasing
            if (is_causal_path && !is_biasing_path) {
                return('causal path')
            } else {
                return('biasing path')
            }
        }
        NA
    }, MARGIN=1)

    attr(nodes, 'dag') = dag
    nodes
}

geom_dagitty_node = function(
    nudge_y=0, legend_label=NULL, implementation='ggstar',
    hide_legend=c(), a=0.15, b=0.1, size=1
) {
    fill_colors = c(
        'exposure'='#d1e04e',
        'outcome'='#4cbde9',
        'ancestor of exposure'='#d1e04e',
        'ancestor of outcome'='#4cbde9',
        'ancestor of exposure and outcome'='#ffabab',
        'adjusted variable'='#fafafa',
        'unobserved (latent)'='#eeeeee',
        'other variable'='#c3c3c3'
    )
    outline_colors = c(
        'exposure'='black',
        'outcome'='black',
        'ancestor of exposure'='#bed403',
        'ancestor of outcome'='#00a2e0',
        'ancestor of exposure and outcome'='#ff7777',
        'adjusted variable'='#050505',
        'unobserved (latent)'='#a7a7a7',
        'other variable'='#666666'
    )
    stopifnot(all(names(fill_colors) == names(outline_colors)))
    legend_limits = setdiff(names(fill_colors), hide_legend)

    if (implementation == 'ggstar') {
        nodes = ggstar::geom_star(
            aes(fill=kind, color=kind, y=y+{{nudge_y}}),
            starshape='ellipse',
            size=size * 5,
            angle=90,
            key_glyph=draw_key_point
        )
    } else {
        nodes = ggforce::geom_ellipse(
            aes(
                fill=kind, color=kind,
                x0=x, y0=y+{{nudge_y}}, a=a*size,
                b=b*size, angle=0
            ),
            key_glyph=draw_key_point
        )
    }

    list(
        nodes,
        geom_point(
            aes(
                shape=ifelse(
                    kind == 'exposure',
                    '▶',
                    ifelse(
                        kind == 'outcome',
                        '▮',
                        NA
                    )
                ),
                y=y+{{nudge_y}}
            ),
            size=2
        ),
        scale_shape_identity(guide='none'),
        guides(
            fill=guide_legend(
                override.aes = list(shape = 'circle filled', size = 5),
            )
        ),
        scale_color_manual(values=outline_colors, limits=legend_limits, name=legend_label),
        scale_fill_manual(values=fill_colors, limits=legend_limits, name=legend_label)
    )
}

geom_node_label = function(nudge_y=0, color='black', bg.colour=alpha('white', 0.8)) {
    shadowtext::geom_shadowtext(
        aes(label=label, y=y+{{nudge_y}}),
        check_overlap = TRUE,
        size=3.5,
        color=color,
        bg.colour=bg.colour
    )
}

daigtty_arrow =  grid::arrow(length = grid::unit(6, "pt"), type = "closed")


draw_line_glyph = function(data, params, size) {
    grid::segmentsGrob(0.1, 0.5, 0.9, 0.5,
        gp = grid::gpar(
            col = alpha(data$edge_colour, data$edge_alpha),
            lwd = 2 * data$edge_width * .pt,
            lty = data$edge_linetype, lineend = 'butt'
        )
    )
}

theme_dagitty = function(background='#eeeeee') {
    list(
        ggdag::theme_dag(),
        theme(panel.background = element_rect(
            fill=background, color=background
        ))
    )
}

get_long_edges = function(gg_edges, dag) {
    function(gg_edges) {
        dag = attr(gg_edges, 'dag')
        dt_edges = dplyr::rename(
            dagitty::edges(dag),
            name=v,
            to=w,
            direction=e,
            xctrl=x,
            yctrl=y
        )
        edges = merge(gg_edges, dt_edges, all=TRUE)
        edges$xctrl = ifelse(
            is.na(edges$xctrl),
            (edges$x + edges$xend) / 2,
            edges$xctrl
        )
        edges$yctrl = ifelse(
            is.na(edges$yctrl),
            (edges$y + edges$yend) / 2,
            edges$yctrl
        )
        edges = reshape(
            edges,
            direction = 'long',
            varying = list(
                c('x', 'xctrl', 'xend'),
                c('y', 'yctrl', 'yend')
            ),
            v.names = c('x', 'y')
        )
        edges
    }
}

node_dagitty_long = function(dag, ...) {
    gg_edges = node_dagitty(dag, ...)
    get_long_edges()(gg_edges)
}

# modified from ggraph::geom_edge_arc2
geom_edge_manual_arc <- function(
    mapping = NULL, data = get_edges('long'),
    position = 'identity', arrow = NULL,
    n = 100, fold = FALSE, lineend = 'butt',
    linejoin = 'round', linemitre = 1,
    label_colour = 'black', label_alpha = 1,
    label_parse = FALSE, check_overlap = FALSE,
    angle_calc = 'rot', force_flip = TRUE,
    label_dodge = NULL, label_push = NULL,
    show.legend = NA, ...
) {
  mapping <-  ggraph:::complete_edge_aes(mapping)
  mapping <- ggraph:::aes_intersect(mapping, aes(
    x = x, y = y, group = edge.id,
    circular = circular
  ))
  layer(
    data = data, mapping = mapping, stat = ggforce::StatBezier2,
    geom = ggraph::GeomEdgePath, position = position, show.legend = show.legend,
    inherit.aes = FALSE,
    params = ggraph:::expand_edge_aes(
      list(
        arrow = arrow, lineend = lineend, linejoin = linejoin,
        linemitre = linemitre, na.rm = FALSE, n = n,
        interpolate = TRUE, fold = fold,
        label_colour = label_colour, label_alpha = label_alpha,
        label_parse = label_parse, check_overlap = check_overlap,
        angle_calc = angle_calc, force_flip = force_flip,
        label_dodge = label_dodge, label_push = label_push, ...
      )
    )
  )
}


# usptream: suggest that both should be arcs but the default curvature for the former should be 0
geom_dag_edges <- function(
    mapping = NULL,
    data_directed = ggdag:::filter_direction("->"),
    data_bidirected = ggdag:::filter_direction("<->"),
    curvature = 0.3,
    directed_curvature=0,
    arrow_directed = grid::arrow(length = grid::unit(5, "pt"), type = "closed"),
    arrow_bidirected = grid::arrow(length = grid::unit(5, "pt"), ends = "both", type = "closed"),
    position = "identity", na.rm = TRUE, show.legend = NA, inherit.aes = TRUE, fold = FALSE,
    ...
) {
    list(
        ggdag::geom_dag_edges_arc(
            mapping,
            data = data_directed,
            arrow = arrow_directed,
            curvature=directed_curvature,
            position = position, na.rm = na.rm,
            show.legend = show.legend, inherit.aes = inherit.aes,
            ...
        ),
        ggdag::geom_dag_edges_arc(
            mapping,
            data = data_bidirected,
            arrow = arrow_bidirected,
            curvature = curvature, position = position,
            na.rm = na.rm, show.legend = show.legend,
            inherit.aes = inherit.aes, fold = fold,
            ...
        )
    )
}

geom_dag_edges2 = function(
    mapping = NULL,
    data = NULL,
    data_directed = ggdag:::filter_direction("->"),
    data_bidirected = ggdag:::filter_direction("<->"),
    arrow_directed = grid::arrow(length = grid::unit(5, "pt"), type = "closed"),
    arrow_bidirected = grid::arrow(length = grid::unit(5, "pt"), ends = "both", type = "closed"),
    start_cap=ggraph::ellipsis(a=1, b=0.5, 'cm'),
    end_cap=start_cap,
    ...
) {
    list(
        geom_edge_manual_arc(
            mapping = mapping,
            data = function (edges) {
                if (is.function(data)) {
                    edges = data(edges)
                } else {
                    edges = data
                }
                if (is.function(data_directed)) {
                    edges = data_directed(edges)
                } else {
                    edges = data_directed
                }
                edges
            },
            arrow = arrow_directed,
            end_cap = end_cap,
            start_cap = start_cap,
            ...
        ),
        geom_edge_manual_arc(
            mapping = mapping,
            data = function (edges) {
                if (is.function(data)) {
                    edges = data(edges)
                } else {
                    edges = data
                }
                if (is.function(data_bidirected)) {
                    edges = data_bidirected(edges)
                } else {
                    edges = data_bidirected
                }
                edges
            },
            arrow = arrow_bidirected,
            end_cap = end_cap,
            start_cap = start_cap,
            ...
        )
    )
}

geom_dagitty_edge = function(legend_label=NULL, ...) {
    colors = c(
        'biasing path'='#d01c8b',
        'causal path'='#4dac26'
    )
    list(
        geom_dag_edges(
            aes(y=y_nudged, yend=yend_nudged, edge_colour=path_kind),
            arrow.fill='white',
            key_glyph=draw_line_glyph,
            ...
        ),
        ggraph::scale_edge_color_manual(
            values=colors, limits=names(colors),
            name=legend_label,
            na.value='black'
        )
    )
}

geom_dagitty_edge2 = function(
    data=get_long_edges(), legend_label=NULL,
    ...
) {
    colors = c(
        'biasing path'='#d01c8b',
        'causal path'='#4dac26'
    )
    list(
        geom_dag_edges2(
            mapping=aes(x=x, group=id, y=y, edge_colour=path_kind),
            data=data,
            arrow.fill='white',
            key_glyph=draw_line_glyph,
            ...
        ),
        ggraph::scale_edge_color_manual(
            values=colors, limits=names(colors),
            name=legend_label,
            na.value='black'
        )
    )
}