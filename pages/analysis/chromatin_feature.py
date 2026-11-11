import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html

from callbacks.constants import Constants

# =====================
# Miscellaneous Modals
# =====================

vcf_format_modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle("Accepted File Format")),
        dbc.ModalBody(
            [
                html.P(
                    "Upload a file containing the variant data for chromatin feature prediction. "
                    "It must include CHR, POS, REF, and ALT columns. "
                ),
                html.P(
                    [
                        "Variant positions should be given with respect to the ",
                        html.B("IRGSP-1.0"),
                        " (Nipponbare) reference genome.",
                    ]
                ),
                html.P(
                    "Only the CHROM, POS, REF, and ALT columns are used; all "
                    "other VCF columns and header metadata are ignored."
                ),
            ]
        ),
    ],
    id="chromatin-file-modal",
    is_open=False,
    size="lg",
    scrollable=True,
)

window_size_modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle("Context Window Size")),
        dbc.ModalBody(
            [
                html.P(
                    "RiSPICE scores each variant using the sequence "
                    "surrounding it. A larger context window gives the "
                    "model more flanking sequence to work with, at the "
                    "cost of diluting the relative influence of the "
                    "mutated position itself."
                ),
                html.P(
                    "Significance thresholds are computed separately for "
                    "each window size against an empirical background "
                    "distribution."
                ),
            ]
        ),
    ],
    id="chromatin-window-modal",
    is_open=False,
    size="lg",
    scrollable=True,
)

# ============
# Main Layout
# ============

layout = html.Div(
    id={"type": "analysis-layout", "label": Constants.LABEL_CHROMATIN},
    hidden=True,
    children=[
        html.Div(
            [
                html.P(
                    [
                        Constants.INTRO_CHROMATIN,
                        " Click ",
                        dcc.Link(
                            [
                                "here ",
                                html.I(
                                    id="demo-link",
                                    className="fa-solid fa-up-right-from-square fa-2xs",
                                ),
                            ],
                            href="https://github.com/bioinfodlsu/rice-pilaf/wiki/2.5-Chromatin-Feature-Prediction",
                            target="_blank",
                            className="top-navbar-item",
                        ),
                        " for the user guide.",
                    ]
                )
            ],
            className="analysis-intro p-3",
        ),
        html.Br(),
        html.Div(
            [
                html.I(className="bi bi-chevron-bar-right me-2 non-clickable"),
                html.Span(id="chromatin-genomic-intervals-input"),
                html.Span("Your Input: Upload a file first.", id="chromatin-genomic-intervals-input"),
            ],
            className="analysis-intro p-3",
        ),
        html.Br(),
        html.Div(
            [
                dbc.Label(
                    [
                        "Upload a file containing the variant data for chromatin feature prediction. It must include CHR, POS, REF, and ALT columns.",
                        html.I(
                            className="bi bi-info-circle",
                            id="chromatin-file-tooltip",
                            n_clicks=0,
                        ),
                    ]
                ),
                vcf_format_modal,
                html.Br(),
                dbc.Alert(
                    id="chromatin-file-error",
                    color="danger",
                    style={"display": "none"},
                ),
                dcc.Upload(
                    id="chromatin-file-upload",
                    children=dbc.Button(
                        [html.I(className="bi bi-cloud-upload me-2"), "Select a file"],
                        color="secondary",
                        outline=True,
                        className="page-button",
                    ),
                    multiple=False,
                    style={"display": "inline-block"},
                ),
                html.Span(
                    "No file selected",
                    id="chromatin-file-filename",
                    className="text-muted ms-3",
                ),
                html.Div(id="chromatin-file-filename", className="mt-2"),
                html.Br(),
                dbc.Label(
                    [
                        "Select the context window size",
                        html.I(
                            className="bi bi-info-circle",
                            id="chromatin-window-tooltip",
                            n_clicks=0,
                        ),
                    ]
                ),
                window_size_modal,
                html.Br(),
                dbc.RadioItems(
                    id="chromatin-window-size",
                    options=[
                        {"value": "500", "label": "500 bp"},
                        {"value": "750", "label": "750 bp"},
                        {"value": "1000", "label": "1 kb"},
                    ],
                    value="1000",
                    inline=True,
                    className="ms-3 mt-1",
                ),
                html.Br(),
                dbc.Button(
                    "Run Analysis",
                    id="chromatin-submit",
                    className="page-button",
                    n_clicks=0,
                ),
            ],
            className="analysis-intro p-3",
        ),
        html.Br(),
        html.Div(
            id="chromatin-results-container",
            style={"display": "none"},
            children=[
                dcc.Loading(
                    id="chromatin-loading",
                    children=[
                        html.Hr(className="mt-3 mb-3"),
                        html.Br(),
                        html.Div(
                            id="chromatin-input", className="analysis-intro p-3"
                        ),
                        html.Br(),
                        html.Div(
                            [
                                html.Div(
                                    [html.Span(id="chromatin-table-stats")],
                                    className="mb-3 stats",
                                ),
                                dbc.Tabs(
                                    id="chromatin-score-tabs",
                                    active_tab="tab-0",
                                    children=[
                                        dcc.Tab(
                                            label="Overall Scores",
                                            value="Overall Scores",
                                        ),
                                        dcc.Tab(
                                            label="Per-Feature Scores",
                                            value="Per-Feature Scores",
                                        ),
                                        dcc.Tab(
                                            label="In Silico Saturated Mutagenesis",
                                            value="In Silico Saturated Mutagenesis",
                                        ),
                                    ],
                                ),
                                html.Br(),
                                html.P(
                                    html.Div(
                                        [
                                            dbc.Button(
                                                [
                                                    html.I(
                                                        className="bi bi-download me-2"
                                                    ),
                                                    "Export to CSV",
                                                ],
                                                id="chromatin-export-table",
                                                n_clicks=0,
                                                color="light",
                                                size="sm",
                                                className="table-button",
                                            ),
                                            dcc.Download(
                                                id="chromatin-download-df-to-csv"
                                            ),
                                            dbc.Button(
                                                [
                                                    html.I(
                                                        className="bi bi-arrow-clockwise me-2"
                                                    ),
                                                    "Reset Table",
                                                ],
                                                id="chromatin-reset-table",
                                                color="light",
                                                size="sm",
                                                className="ms-3 table-button",
                                            ),
                                        ],
                                        style={"textAlign": "right"},
                                    )
                                ),
                                dash_table.DataTable(
                                    id="chromatin-scores-table",
                                    style_cell={"whiteSpace": "pre-line"},
                                    markdown_options={"html": True},
                                    sort_action="native",
                                    filter_action="native",
                                    filter_options={
                                        "case": "insensitive",
                                        "placeholder_text": "🔎︎ Search Column",
                                    },
                                    page_action="native",
                                    page_size=15,
                                    cell_selectable=False,
                                    style_table={"overflowX": "auto"},
                                ),
                                html.Br(),
                                html.Br(),
                                html.Div(
                                    id="chromatin-mutagenesis-container",
                                    children=[
                                        html.P(
                                            "The heatmap below shows the predicted "
                                            "regulatory effect of every possible "
                                            "single-nucleotide substitution within "
                                            "the selected window, centered on your "
                                            "variant of interest."
                                        ),
                                        html.Div(id="chromatin-mutagenesis-heatmap"),
                                    ],
                                    className="analysis-intro p-3",
                                ),
                            ],
                            className="analysis-intro p-3",
                        ),
                    ],
                )
            ],
        ),
    ],
    className="mt-2 mb-4",
)
