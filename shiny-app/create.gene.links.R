#'
#'
create.group.containers <- function (genes) {
    gene.links.container.id <- "gene-links"
    gene.containers <- mapply(function (gene.group, name, seq) {
        group.id <- sprintf("gr-%s-%07d", gsub("[^-0-9a-zA-Z]+", "-", name), sample.int(9999999L, 1L))

        group.sep <- div(
            class = "panel-heading",
            h4(
                class = "panel-title",
                actionLink(
                    sprintf("gene-group-%d", seq),
                    label = sprintf("%s", name),
                    icon = icon("arrow-circle-right", class = "right"),
                    href = sprintf("#%s", group.id),
                    `data-toggle` = "collapse",
                    `data-parent` = sprintf("#%s", gene.links.container.id)
                )
            )
        )

        list.group.collapsible <- div(
            id = group.id,
            class = "panel-collapse collapse",
            loaded = FALSE,
            collapsed = TRUE
        )

        div(
            class = "panel panel-default",
            group.sep,
            list.group.collapsible
        )
    }, genes, names(genes), seq_along(genes), SIMPLIFY = FALSE)

    containers <- do.call(div, c(unname(gene.containers),
                                 list(class = "panel-group", id = gene.links.container.id)))

    gene.links <- mapply(function (gene.group, name, gr.seq) {
        links <- mapply(function (gene, gene.seq) {
            actionLink(
                sprintf("gene-id-%d-%d", gr.seq, gene.seq),
                label = gene,
                class = "list-group-item"
            )
        }, gene.group, seq_along(gene.group), SIMPLIFY = FALSE)

        do.call(div, c(unname(links), list(class = "list-group")))
    }, genes, names(genes), seq_along(genes), SIMPLIFY = FALSE)


    return(list(
        containers = containers,
        gene.links = gene.links
    ))
}
