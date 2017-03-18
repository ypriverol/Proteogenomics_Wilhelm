Shiny.addCustomMessageHandler("scrollCallback", function (msg) {
    window.scrollTo(0, 0);
});


function renderGene (item, escape) {
    return '<div class="gene-item">' +
        '<span class="gene-description">' +
            '<span class="gene-symbol">' + escape(item.gene) + '</span>' +
            '<span class="gene-ensembl">' + escape(item.ensembl) + '</span>' +
        '</span>' +
        '<span class="gene-name">' + escape(item.geneName) + '</span>' +
        '<ul class="gene-stats">' +
            '<li>Cor = ' + item.cor + '</li>' +
            '<li># Observations: ' + item.availBoth + '</li>' +
        '</ul">' +
    '</div>';
}
