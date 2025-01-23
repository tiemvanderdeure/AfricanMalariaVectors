using SummaryTables, DataFrames, Printf
import WriteDocx as W

### Table with evaluation scores
get_mean_score(e::SDM.SDMensembleEvaluation, t, m) = mean(e.results[t][m].score)
scores = [
    get_mean_score(e, t, m)
    for e in evaluations, m in keys(measures), t in (:train, :test)
]

species_cells = Cell.("A. " .* as_label.(FOCUS_SPECIES), italic = true, halign = :left) |> collect
scores_cells = mapreduce(hcat, keys(measures)) do m
    map(collect(evaluations)) do e
        Cell(@sprintf "%0.2f (%0.2f)" get_mean_score(e, :test, m) get_mean_score(e, :train, m))
    end
end
col_titles = Cell.(["Species"; collect(string.(keys(measures)))], bold = true, border_bottom = true)

cells = [col_titles'; [species_cells scores_cells]]
table1 = SummaryTables.Table(cells)

### Table with variable importances
varimp_cells = mapreduce(vcat, collect(variable_importances)) do varimps
    [Cell(@sprintf "%0.3f" i) for i in varimps]'
end
colnames = Cell.(["Species"; collect(string.(keys(first(variable_importances))))], bold = true, border_bottom = true)
cells = [colnames'; [species_cells varimp_cells]]
table2 = SummaryTables.Table(cells)


doc = W.Document(
    W.Body([
        W.Section([
            W.Paragraph([
                W.Run([W.Text("Table 1")]),
            ]),
            SummaryTables.to_docx(table1),
            W.Paragraph([
                W.Run([W.Text("Table 2")]),
            ]),
            SummaryTables.to_docx(table2),
        ]),
    ]),
)

W.save(joinpath("tables", "tables_natcc.docx"), doc)