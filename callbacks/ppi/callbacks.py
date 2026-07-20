import pickle
from collections import namedtuple

import pandas as pd
from dash import Input, Output, State, ctx, dcc, html
from dash.exceptions import PreventUpdate

from ..constants import Constants
from ..general_util import NULL_PLACEHOLDER, get_num_unique_entries, purge_html_export_table
from ..homepage.util import clear_specific_dccStore_data
from ..lift_over.util import get_genes_in_Nb, get_genomic_intervals_from_input, is_error
from ..links_util import (
    get_msu_browser_link_single_str,
    get_rgi_orthogroup_link_single_str,
)
# The module-detection/enrichment/graph-display engine is shared across the
# Coexpression and PPI modules (both operate over gene networks precomputed
# by the same offline pipeline); only the network catalog and the
# PPI-specific topological statistics/hub-gene logic are PPI-specific
from ..coexpression.util import (
    check_if_valid_msu_ids,
    convert_to_df,
    count_genes_in_module,
    count_modules,
    do_module_enrichment_analysis,
    get_gene_description_entry,
    get_interpro_entry,
    get_noun_for_active_tab,
    get_parameters_for_algo,
    get_pfam_entry,
    get_pubmed_entry,
    get_qtaro_entry,
    get_rapdb_entry,
    get_uniprot_entry,
    get_user_facing_algo,
    get_user_facing_parameter,
    load_module_graph,
    module_detection_algos,
    sanitize_msu_id,
)
from .util import get_query_network_summary, get_user_facing_network

Parameter_slider = namedtuple("Parameter_slider", ["marks", "value"])


def init_callback(app):
    @app.callback(
        Output("ppi-genomic-intervals-input", "children"),
        State("homepage-submitted-genomic-intervals", "data"),
        Input("homepage-is-submitted", "data"),
        Input("ppi-submit", "n_clicks"),
    )
    def display_input(nb_intervals_str, homepage_is_submitted, *_):
        """
        Displays the genomic interval input in the ppi page

        Parameters:
        - nb_intervals_str: Submitted genomic interval
        - homepage_is_submitted: [Homepage] Saved boolean value of submitted valid input
        - *_: Other input that facilitates displaying of the submitted genomic interval

        Returns:
        - Submitted genomic interval text
        """

        if homepage_is_submitted:
            if nb_intervals_str and not is_error(
                get_genomic_intervals_from_input(nb_intervals_str)
            ):
                return [html.B("Your Input Intervals: "), html.Span(nb_intervals_str)]

            return None

        raise PreventUpdate

    # =================
    # Input-related
    # =================
    @app.callback(
        Output("session-container", "children", allow_duplicate=True),
        Output("ppi-is-submitted", "data", allow_duplicate=True),
        Output("ppi-submitted-addl-genes", "data", allow_duplicate=True),
        Output("ppi-valid-addl-genes", "data", allow_duplicate=True),
        Output("ppi-combined-genes", "data", allow_duplicate=True),
        Output("ppi-submitted-network", "data", allow_duplicate=True),
        Output("ppi-submitted-clustering-algo", "data", allow_duplicate=True),
        Output("ppi-submitted-parameter-slider", "data", allow_duplicate=True),
        Output("ppi-addl-genes-error", "style", allow_duplicate=True),
        Output("ppi-addl-genes-error", "children", allow_duplicate=True),
        Input("ppi-submit", "n_clicks"),
        State("homepage-is-submitted", "data"),
        State("session-container", "children"),
        State("homepage-submitted-genomic-intervals", "data"),
        State("ppi-addl-genes", "value"),
        State("ppi-network", "value"),
        State("ppi-clustering-algo", "value"),
        State("ppi-parameter-slider", "marks"),
        State("ppi-parameter-slider", "value"),
        prevent_initial_call=True,
    )
    def submit_ppi_input(
        ppi_submit_n_clicks,
        homepage_is_submitted,
        dccStore_children,
        genomic_intervals,
        submitted_addl_genes,
        submitted_network,
        submitted_algo,
        submitted_slider_marks,
        submitted_slider_value,
    ):
        """
        Parses ppi input, displays the ppi result container
        - If user clicks on the ppi submit button, the inputs will be parsed and either an error message or the ppi results container will appear

        Parameters:
        - ppi_submit_n_clicks: Number of clicks pressed on the tfbs submit button
        - homepage_is_submitted: [Homepage] Saved boolean value of submitted valid input
        - dccStore_children: List of dcc.Store data
        - genomic_intervals: Saved genomic intervals found in the dcc.Store
        - submitted_addl_genes: Submitted ppi additional genes
        - submitted_network: Submitted ppi network
        - submitted_algo: Submitted ppi clustering algorithm
        - submitted_slider_marks: Submitted parameter slider marks
        - submitted_slider_value: Submitted parameter slider value

        Returns:
        - ('session-container', 'children'): Updated dcc.Store data
        - ('ppi-is-submitted', 'data'): [PPI] True for submitted valid input; otherwise False
        - ('ppi-submitted-addl-genes', 'data'): Submitted ppi additional genes
        - ('ppi-valid-addl-genes', 'data'): Submitted ppi valid additional genes
        - ('ppi-combined-genes', 'data'): PPI combined genes
        - ('ppi-submitted-network', 'data'): Submitted ppi network
        - ('ppi-submitted-clustering-algo', 'data'): Submitted ppi clustering algorithm
        - ('ppi-submitted-parameter-slider', 'data): Submitted ppi parameter slider tuple
        - ('ppi-addl-genes-error', 'style'): {'display': 'block'} for displaying the error message; otherwise {'display': 'none'}
        - ('ppi-addl-genes-error', 'children'): Error message
        """

        if homepage_is_submitted and ppi_submit_n_clicks >= 1:
            parameter_slider_value = Parameter_slider(
                submitted_slider_marks, submitted_slider_value
            )._asdict()
            submitted_parameter_slider = {submitted_algo: parameter_slider_value}

            if submitted_addl_genes:
                submitted_addl_genes = submitted_addl_genes.strip()
            else:
                submitted_addl_genes = ""

            list_addl_genes = list(
                filter(
                    None,
                    [
                        sanitize_msu_id(gene.strip())
                        for gene in submitted_addl_genes.split(";")
                    ],
                )
            )

            # Check which genes are valid MSU IDs
            list_addl_genes, invalid_genes = check_if_valid_msu_ids(list_addl_genes)

            if not invalid_genes:
                error_display = {"display": "none"}
                error = None
            else:
                error_display = {"display": "block"}

                if len(invalid_genes) == 1:
                    error_msg = invalid_genes[0] + " is not a valid MSU accession ID."
                    error_msg_ignore = "It"
                else:
                    if len(invalid_genes) == 2:
                        error_msg = invalid_genes[0] + " and " + invalid_genes[1]
                    else:
                        error_msg = (
                            ", ".join(invalid_genes[:-1]) + ", and " + invalid_genes[-1]
                        )

                    error_msg += " are not valid MSU accession IDs."
                    error_msg_ignore = "They"

                error = [
                    html.Span(error_msg),
                    html.Br(),
                    html.Span(
                        f"{error_msg_ignore} will be ignored when running the analysis."
                    ),
                ]

            # Perform lift-over if it has not been performed.
            # Otherwise, just fetch the results from the file
            implicated_gene_ids = get_genes_in_Nb(genomic_intervals)[1]

            gene_ids = list(set.union(set(implicated_gene_ids), set(list_addl_genes)))

            dccStore_children = clear_specific_dccStore_data(
                dccStore_children,
                "ppi-pathway-active",
                "ppi-graph-active",
            )

            return (
                dccStore_children,
                True,
                submitted_addl_genes,
                list_addl_genes,
                gene_ids,
                submitted_network,
                submitted_algo,
                submitted_parameter_slider,
                error_display,
                error,
            )

        raise PreventUpdate

    @app.callback(
        Output("ppi-results-container", "style"),
        Input("ppi-is-submitted", "data"),
    )
    def display_ppi_output(ppi_is_submitted):
        """
        Displays the ppi results container

        Parameters:
        - ppi_is_submitted: [PPI] Saved boolean value of submitted valid input

        Returns:
        - ('ppi-results-container', 'style'): {'display': 'block'} for displaying the ppi results container; otherwise {'display': 'none'}
        """

        if ppi_is_submitted:
            return {"display": "block"}

        else:
            return {"display": "none"}

    @app.callback(
        Output("ppi-submit", "disabled"),
        Input("ppi-submit", "n_clicks"),
        Input("ppi-module-graph", "elements"),
        Input("ppi-pathways", "data"),
        Input("ppi-module-stats", "children"),
    )
    def disable_ppi_button_upon_run(n_clicks, *_):
        """
        Disables the submit button in the ppi page until computation is done in the ppi page

        Parameters:
        - n_clicks: Number of clicks pressed on the ppi submit button
        - *_: Other input that facilitates the disabling of the coexpressino submit button

        Returns:
        - ('ppi-submit', 'disabled'): True for disabling the submit button; otherwise False
        """

        return ctx.triggered_id == "ppi-submit" and n_clicks > 0

    @app.callback(
        Output("ppi-clustering-algo-modal", "is_open"),
        Output("ppi-network-modal", "is_open"),
        Output("ppi-parameter-modal", "is_open"),
        Output("ppi-converter-modal", "is_open"),
        Input("ppi-clustering-algo-tooltip", "n_clicks"),
        Input("ppi-network-tooltip", "n_clicks"),
        Input("ppi-parameter-tooltip", "n_clicks"),
        Input("ppi-converter-tooltip", "n_clicks"),
    )
    def open_modals(
        algo_tooltip_n_clicks,
        network_tooltip_n_clicks,
        parameter_tooltip_n_clicks,
        converter_tooltip_n_clicks,
    ):
        """
        Displays the ppi tooltip modals

        Parameters:
        - algo_tooltip_n_clicks: Number of clicks pressed for the tooltip button near the ppi clustering algorithm input field
        - network_tooltip_n_clicks: Number of clicks pressed for the tooltip button near the ppi network input field
        - parameter_tooltip_n_clicks: Number of clicks pressed for the tooltip button near the ppi parameter slider input field
        - converter_tooltip_n_clicks:Number of clicks pressed for the tooltip button near the ppi additional genes input field

        Returns:
        - ('ppi-clustering-algo-modal', 'is_open'): True for showing the ppi clustering algorithm tooltip; otherwise False
        - ('ppi-network-modal', 'is_open'): True for showing the ppi network tooltip; otherwise False
        - ('ppi-parameter-modal', 'is_open'): True for showing the ppi parameter slider tooltip; otherwise False
        - ('ppi-converter-modal', 'is_open'): True for showing the tfbs additional genes tooltip; otherwise False
        """

        if (
            ctx.triggered_id == "ppi-clustering-algo-tooltip"
            and algo_tooltip_n_clicks > 0
        ):
            return True, False, False, False

        if (
            ctx.triggered_id == "ppi-network-tooltip"
            and network_tooltip_n_clicks > 0
        ):
            return False, True, False, False

        if (
            ctx.triggered_id == "ppi-parameter-tooltip"
            and parameter_tooltip_n_clicks > 0
        ):
            return False, False, True, False

        if (
            ctx.triggered_id == "ppi-converter-tooltip"
            and converter_tooltip_n_clicks > 0
        ):
            return False, False, False, True

        raise PreventUpdate

    @app.callback(
        Output("ppi-parameter-slider", "marks"),
        Output("ppi-parameter-slider", "value"),
        Input("ppi-clustering-algo", "value"),
        State("ppi-submitted-parameter-slider", "data"),
    )
    def set_parameter_slider(algo, parameter_slider):
        """
        Sets the parameter slider data

        Parameters:
        - algo: Selected clustering algorithm
        - parameter_slider: Selected value of the parameter slider

        Returns:
        - ('ppi-parameter-slider', 'marks'): Parameter slider marks
        - ('ppi-parameter-slider', 'value'): Parameter slider value
        """

        if parameter_slider and algo in parameter_slider:
            return parameter_slider[algo]["marks"], parameter_slider[algo]["value"]

        try:
            return (
                get_parameters_for_algo(algo, network="STRING-Physical"),
                module_detection_algos[algo].default_param
                * module_detection_algos[algo].multiplier,
            )
        except FileNotFoundError:
            # STOPGAP: the module-detection pipeline hasn't been run locally
            # for this algorithm/network yet (see Constants.NETWORK_MODULES).
            # Degrade gracefully instead of crashing the callback so the rest
            # of the page (inputs, network stats, hub genes) stays usable.
            return {0: "Data not yet available"}, 0

    @app.callback(
        Output("ppi-input", "children"),
        Input("ppi-is-submitted", "data"),
        State("ppi-valid-addl-genes", "data"),
        State("ppi-submitted-network", "data"),
        State("ppi-submitted-clustering-algo", "data"),
        State("ppi-submitted-parameter-slider", "data"),
    )
    def display_ppi_submitted_input(
        ppi_is_submitted, genes, network, algo, submitted_parameter_slider
    ):
        """
        Displays the ppi submitted input

        Parameters:
        - ppi_is_submitted: [PPI] Saved boolean value of submitted valid input
        - genes: Saved ppi valid additional genes found in the dcc.Store
        - network: Saved ppi network found in the dcc.Store
        - submitted_parameter_slider: Saved ppi parameter slider tuple found in the dcc.Store

        Returns:
        - ('ppi-input', 'children'): Submitted ppi inputs text
        """

        if ppi_is_submitted:
            parameters = 0
            if submitted_parameter_slider and algo in submitted_parameter_slider:
                parameters = submitted_parameter_slider[algo]["value"]

            if not genes:
                genes = "None"
            else:
                genes = "; ".join(genes)

            return [
                html.B("Additional Genes: "),
                genes,
                html.Br(),
                html.B("Selected PPI Network: "),
                get_user_facing_network(network),
                html.Br(),
                html.B("Selected Module Detection Algorithm: "),
                get_user_facing_algo(algo),
                html.Br(),
                html.B("Selected Algorithm Parameter: "),
                get_user_facing_parameter(algo, parameters, network="STRING-Physical"),
            ]

        raise PreventUpdate

    # =================
    # Module-related
    # =================

    @app.callback(
        Output("ppi-modules", "options"),
        Output("ppi-modules", "value"),
        Output("ppi-results-module-tabs-container", "style"),
        Output("ppi-module-stats", "children"),
        State("homepage-submitted-genomic-intervals", "data"),
        Input("ppi-combined-genes", "data"),
        State("ppi-valid-addl-genes", "data"),
        Input("ppi-submitted-network", "data"),
        Input("ppi-submitted-clustering-algo", "data"),
        State("homepage-is-submitted", "data"),
        State("ppi-submitted-parameter-slider", "data"),
        State("ppi-submitted-module", "data"),
        State("ppi-is-submitted", "data"),
    )
    def perform_module_enrichment(
        genomic_intervals,
        combined_gene_ids,
        valid_addl_genes,
        submitted_network,
        submitted_algo,
        homepage_is_submitted,
        submitted_parameter_slider,
        module,
        ppi_is_submitted,
    ):
        """
        Displays the ppi pathways table

        Parameters:
        - genomic_intervals: Saved genomic intervals found in the dcc.Store
        - combined_gene_ids: Saved combined gene ids found in the dcc.Store
        - valid_addl_genes: Saved ppi valid additional genes found in the dcc.Store
        - submitted_network: Saved ppi network found in the dcc.Store
        - submitted_algo: Saved ppi clustering algorithm found in the dcc.Store
        - homepage_is_submitted: [Homepage] Saved boolean value of submitted valid input
        - submitted_parameter_slider: Saved parameters slider tuple found in the dcc.Store
        - module: Saved selected ppi module found in the dcc.Store
        - ppi_is_submitted: [PPI] Saved boolean value of submitted valid input

        Returns:
        - ('ppi-modules', 'options'): List of available modules
        - ('ppi-modules', 'value'): Selected module value
        - ('ppi-results-module-tabs-container', 'style'): {'display': 'block'} for displaying the module tabs container; otherwise {'display': 'none'}
        - ('ppi-module-stats', 'children'): Stats for the ppi module
        """

        if homepage_is_submitted and ppi_is_submitted:
            if submitted_algo and submitted_algo in submitted_parameter_slider:
                parameters = submitted_parameter_slider[submitted_algo]["value"]

                try:
                    enriched_modules = do_module_enrichment_analysis(
                        combined_gene_ids,
                        genomic_intervals,
                        valid_addl_genes,
                        submitted_network,
                        submitted_algo,
                        parameters,
                    )
                    total_num_modules = count_modules(
                        submitted_network, submitted_algo, parameters
                    )
                except FileNotFoundError as e:
                    # STOPGAP: no module-detection output exists locally yet
                    # for this network/algorithm (pipeline hasn't been run,
                    # or its results haven't been synced from the lab
                    # server). Show a clear message instead of crashing.
                    # TEMP DEBUG: surface the exact missing path so we can
                    # pin down the real cause instead of guessing further.
                    stats = (
                        "No module data available locally for "
                        f"{get_user_facing_network(submitted_network)} / "
                        f"{get_user_facing_algo(submitted_algo)} yet. "
                        "Run the data-prep pipeline for this network/algorithm, "
                        "or sync the generated data from the lab server. "
                        f"[DEBUG: missing file was {e.filename}, "
                        f"parameters={parameters!r} (type={type(parameters).__name__})]"
                    )
                    return [], None, {"display": "none"}, stats

                # Display statistics
                num_enriched_modules = len(enriched_modules)
                stats = f"{num_enriched_modules} out of {total_num_modules} "
                if total_num_modules == 1:
                    stats += "module "
                else:
                    stats += "modules "

                if num_enriched_modules == 1:
                    stats += "was "
                else:
                    stats += "were "

                stats += "found to be enriched (adjusted p-value < 0.05)."

                first_module = None
                if enriched_modules:
                    first_module = enriched_modules[0]
                    module = first_module
                else:
                    return enriched_modules, first_module, {"display": "none"}, stats

                if module:
                    first_module = module

                return enriched_modules, first_module, {"display": "block"}, stats

        raise PreventUpdate

    # =================
    # Table-related
    # =================

    @app.callback(
        Output("ppi-pathways", "data"),
        Output("ppi-pathways", "columns"),
        Output("ppi-graph-stats", "children"),
        Output("ppi-table-stats", "children"),
        Output("ppi-table-container", "style"),
        Input("ppi-combined-genes", "data"),
        Input("ppi-submitted-network", "data"),
        Input("ppi-submitted-clustering-algo", "data"),
        Input("ppi-modules-pathway", "active_tab"),
        Input("ppi-modules", "value"),
        State("ppi-submitted-parameter-slider", "data"),
        State("ppi-is-submitted", "data"),
    )
    def display_pathways(
        combined_gene_ids,
        submitted_network,
        submitted_algo,
        active_tab,
        module,
        submitted_parameter_slider,
        ppi_is_submitted,
    ):
        """
        Displays the ppi pathways table

        Parameters:
        - combined_gene_ids: Saved combined gene ids found in the dcc.Store
        - submitted_network: Saved ppi network found in the dcc.Store
        - submitted_algo: Saved ppi clustering algorithm found in the dcc.Store
        - active_tab: Active tab for a specific ppi table
        - module: Selected ppi module
        - submitted_parameter_slider: Saved parameters slider tuple found in the dcc.Store
        - ppi_is_submitted: [PPI] Saved boolean value of submitted valid input

        Returns:
        - ('ppi-pathways', 'data'): Data for the ppi table depending on the active tab
        - ('ppi-pathways', 'columns'): List of columns for a specific ppi table
        - ('ppi-graph-stats', 'children'): Stats for the ppi graph
        - ('ppi-table-stats', 'children'): Stats for the ppi table
        - ('ppi-table-container', 'style'): {'visibility': 'visible'} for displaying the table container; otherwise {'display': 'none'}
        """

        if ppi_is_submitted:
            if (
                submitted_network
                and submitted_algo
                and submitted_algo in submitted_parameter_slider
            ):
                parameters = submitted_parameter_slider[submitted_algo]["value"]

                try:
                    module_idx = module.split(" ")[1]
                    table, _ = convert_to_df(
                        active_tab,
                        module_idx,
                        submitted_network,
                        submitted_algo,
                        parameters,
                    )
                except Exception:
                    table, _ = convert_to_df(
                        active_tab, None, submitted_network, submitted_algo, parameters
                    )

                columns = [
                    {"id": x, "name": x, "presentation": "markdown"}
                    for x in table.columns
                ]

                num_enriched = get_num_unique_entries(table, "ID")
                stats = f"This module is enriched in {num_enriched} "
                if num_enriched == 1:
                    stats += f"{get_noun_for_active_tab(active_tab).singular}."
                else:
                    stats += f"{get_noun_for_active_tab(active_tab).plural}."

                graph_stats = "The selected module has "
                try:
                    total_num_genes, num_combined_gene_ids = count_genes_in_module(
                        combined_gene_ids,
                        int(module_idx),
                        submitted_network,
                        submitted_algo,
                        parameters,
                    )
                except UnboundLocalError:
                    total_num_genes, num_combined_gene_ids = 0, 0

                if total_num_genes == 1:
                    graph_stats += f"{total_num_genes} gene"
                else:
                    graph_stats += f"{total_num_genes} genes"

                graph_stats += f", among which {num_combined_gene_ids} "

                if num_combined_gene_ids == 1:
                    graph_stats += "is "
                else:
                    graph_stats += "are "
                graph_stats += "implicated by your GWAS/QTL or among those that you manually added."

                if total_num_genes == 0:
                    return (
                        table.to_dict("records"),
                        columns,
                        graph_stats,
                        stats,
                        {"display": "none"},
                    )
                else:
                    return (
                        table.to_dict("records"),
                        columns,
                        graph_stats,
                        stats,
                        {"visibility": "visible"},
                    )

        raise PreventUpdate

    @app.callback(
        Output("ppi-pathways", "filter_query"),
        Output("ppi-pathways", "page_current"),
        Input("ppi-reset-table", "n_clicks"),
        Input("ppi-submit", "n_clicks"),
        Input("ppi-modules-pathway", "active_tab"),
        Input("ppi-modules", "value"),
    )
    def reset_table_filter_page(*_):
        """
        Resets the ppi table and the current page to its original state

        Parameters:
        - *_: Other input that facilitates the resetting of the ppi table

        Returns:
        - ('ppi-pathways', 'filter_query'): '' for removing the filter query
        - ('ppi-pathways', 'page_current'): 0
        """

        return "", 0

    @app.callback(
        Output("ppi-table-container", "style", allow_duplicate=True),
        Input("ppi-submit", "n_clicks"),
        prevent_initial_call=True,
    )
    def hide_table(*_):
        """
        Hides the ppi table

        Parameters:
        - *_: Other inputs to facilitate the hiding of the ppi table

        Returns:
        - ('ppi-table-container', 'style'): {'visibility': 'hidden'} for hiding the ppi table
        """

        return {"visibility": "hidden"}

    @app.callback(
        Output("ppi-download-df-to-csv", "data"),
        Input("ppi-export-table", "n_clicks"),
        State("ppi-pathways", "data"),
        State("ppi-modules", "value"),
    )
    def download_ppi_table_to_csv(download_n_clicks, ppi_df, module):
        """
        Export the ppi table in csv file format

        Parameters:
        - download_n_clicks: Number of clicks pressed on the export ppi table button
        - ppi_df: ppi table data in dataframe format
        - module: Selected ppi module

        Returns:
        - ('ppi-download-df-to-csv', 'data'): PPI table in csv file format data
        """

        if download_n_clicks >= 1:
            df = pd.DataFrame(purge_html_export_table(ppi_df))
            return dcc.send_data_frame(
                df.to_csv,
                f"[{module}] PPI Network Analysis Table.csv",
                index=False,
            )

        raise PreventUpdate

    # =================
    # Graph-related
    # =================

    @app.callback(
        Output("ppi-module-graph", "elements", allow_duplicate=True),
        Output("ppi-module-graph", "layout", allow_duplicate=True),
        Output("ppi-module-graph", "style", allow_duplicate=True),
        Output("ppi-graph-container", "style", allow_duplicate=True),
        Output("ppi-module-graph-node-data", "children", allow_duplicate=True),
        Output(
            "ppi-module-graph-node-data-container",
            "style",
            allow_duplicate=True,
        ),
        Input("ppi-combined-genes", "data"),
        Input("ppi-modules", "value"),
        State("ppi-submitted-network", "data"),
        State("ppi-submitted-clustering-algo", "data"),
        State("ppi-submitted-parameter-slider", "data"),
        Input("ppi-graph-layout", "value"),
        State("ppi-is-submitted", "data"),
        State("ppi-modules", "options"),
        Input("ppi-reset-graph", "n_clicks"),
        prevent_initial_call=True,
    )
    def display_graph(
        combined_gene_ids,
        module,
        submitted_network,
        submitted_algo,
        submitted_parameter_slider,
        layout,
        ppi_is_submitted,
        modules,
        *_,
    ):
        """
        Displays the ppi graph

        Parameters:
        - combined_gene_ids: Saved ppi combined genes found in the dcc.Store
        - module: Selected ppi module
        - submitted_network: Saved ppi network found in the dcc.Store
        - submitted_algo: Saved ppi clustering algorithm in the dcc.Store
        - submitted_parameter_slider: Saved parameter slider tuple found in the dcc.Store
        - layout: Selected ppi graph layout
        - ppi_is_submitted: [PPI] Saved boolean value of submitted valid input
        - modules: List of available modules
        - *_: Other inputs that facilitates the state of the ppi graph

        Returns:
        - ('ppi-module-graph', 'elements'): List of elements of the ppi graph
        - ('ppi-module-graph', 'layout'): Selected ppi graph layout
        - ('ppi-module-graph', 'style'): {'visibility': 'visible'} for displaying the ppi graph; otherwise {'display': 'none'}
        - ('ppi-graph-container', 'style'): {'visibility': 'visible'} for displaying the ppi graph container; otherwise {'display': 'none'}
        - ('ppi-module-graph-node-data', 'children'): Short instruction on how to display the selected node data
        - ('ppi-module-graph-node-data-container', 'style'): {'display': 'block'} for displaying the selected node data; otherwise {'display': 'none'}
        """

        if ppi_is_submitted:
            if (
                submitted_network
                and submitted_algo
                and submitted_algo in submitted_parameter_slider
            ):
                parameters = submitted_parameter_slider[submitted_algo]["value"]

                if not modules:
                    module_graph = load_module_graph(
                        combined_gene_ids,
                        "Click on a node to display information about the gene.",
                        submitted_network,
                        submitted_algo,
                        parameters,
                        layout,
                    )
                else:
                    module_graph = load_module_graph(
                        combined_gene_ids,
                        module,
                        submitted_network,
                        submitted_algo,
                        parameters,
                        layout,
                    )

                # No enriched modules
                if not modules:
                    return module_graph + ({"display": "none"}, "", {"display": "none"})

                return module_graph + (
                    {"visibility": "visible", "width": "100%"},
                    "Click on a node to display information about the gene.",
                    {"display": "block"},
                )

        raise PreventUpdate

    @app.callback(
        Output("ppi-module-graph-node-data", "children"),
        Input("ppi-module-graph", "tapNodeData"),
    )
    def display_node_data(node_data):
        """
        Displays the selected ppi graph's node's data

        Parameters:
        - node_data: Selected ppi graph's node's data

        Returns:
        - ('ppi-module-graph-node-data', 'children'): Selected node data
        """

        if node_data:
            with open(
                f"{Constants.OGI_MAPPING}/Nb_to_ogi.pickle", "rb"
            ) as ogi_file, open(Constants.QTARO_DICTIONARY, "rb") as qtaro_file, open(
                f"{Constants.IRIC}/interpro.pickle", "rb"
            ) as interpro_file, open(
                f"{Constants.IRIC}/pfam.pickle", "rb"
            ) as pfam_file, open(
                f"{Constants.IRIC_MAPPING}/msu_to_iric.pickle", "rb"
            ) as iric_mapping_file, open(
                f"{Constants.TEXT_MINING_PUBMED}", "rb"
            ) as pubmed_file, open(
                f"{Constants.MSU_MAPPING}/msu_to_rap.pickle", "rb"
            ) as rapdb_file, open(
                f"{Constants.GENE_DESCRIPTIONS}/Nb/Nb_gene_descriptions.pickle", "rb"
            ) as gene_descriptions_file:
                ogi_mapping = pickle.load(ogi_file)
                qtaro_mapping = pickle.load(qtaro_file)
                interpro_mapping = pickle.load(interpro_file)
                pfam_mapping = pickle.load(pfam_file)
                iric_mapping = pickle.load(iric_mapping_file)
                pubmed_mapping = pickle.load(pubmed_file)
                rapdb_mapping = pickle.load(rapdb_file)
                gene_descriptions_mapping = pickle.load(gene_descriptions_file)

                gene = node_data["id"]

                node_data = [
                    html.H5("Gene Information", className="pb-3"),
                    html.B("Name: "),
                    get_msu_browser_link_single_str(gene, dash=True),
                    html.Br(),
                    html.B("OGI: "),
                    get_rgi_orthogroup_link_single_str(ogi_mapping[gene], dash=True),
                    html.Br(),
                    html.B(
                        "RAP-DB: ",
                    ),
                    get_rapdb_entry(gene, rapdb_mapping),
                    html.Br(),
                    html.B("Description: "),
                    get_gene_description_entry(gene, gene_descriptions_mapping),
                    html.Br(),
                    html.B("UniProtKB/Swiss-Prot: "),
                    get_uniprot_entry(gene, gene_descriptions_mapping),
                    html.Br(),
                    html.Br(),
                    html.B("Pfam: "),
                    get_pfam_entry(gene, pfam_mapping, iric_mapping),
                    html.Br(),
                    html.B("InterPro: "),
                    get_interpro_entry(gene, interpro_mapping, iric_mapping),
                    html.Br(),
                    html.B("QTL Analyses: "),
                    get_qtaro_entry(gene, qtaro_mapping),
                    html.Br(),
                    html.B("PubMed Article IDs: "),
                    get_pubmed_entry(gene, pubmed_mapping),
                ]

                return node_data

        raise PreventUpdate

    @app.callback(
        Output("ppi-module-graph", "style", allow_duplicate=True),
        Input("ppi-modules", "value"),
        prevent_initial_call=True,
    )
    def hide_graph(*_):
        """
        Hides the ppi graph

        Parameters:
        - *_: Other inputs to facilitate the hiding of the ppi graph

        Returns:
        - ('ppi-module-graph', 'style'): {'visibility': 'hidden'} for hiding the ppi graph
        """

        return {"visibility": "hidden"}

    @app.callback(
        Output("ppi-module-graph", "elements"),
        Output("ppi-module-graph", "layout"),
        Output("ppi-module-graph", "style", allow_duplicate=True),
        Output("ppi-graph-container", "style"),
        Input("ppi-combined-genes", "data"),
        Input("ppi-submitted-network", "data"),
        Input("ppi-submitted-clustering-algo", "data"),
        State("ppi-is-submitted", "data"),
        State("ppi-submitted-parameter-slider", "data"),
        State("ppi-graph-active-layout", "data"),
        prevent_initial_call=True,
    )
    def hide_table_graph(
        combined_gene_ids,
        submitted_network,
        submitted_algo,
        ppi_is_submitted,
        submitted_parameter_slider,
        layout,
    ):
        """
        Hides the ppi graph

        Parameters:
        - combined_gene_ids: Saved ppi combined genes found in the dcc.Store
        - submitted_network: Saved ppi network found in the dcc.Store
        - submitted_algo: Saved ppi clustering algorithm in the dcc.Store
        - ppi_is_submitted: [PPI] Saved boolean value of submitted valid input
        - submitted_parameter_slider: Saved parameter slider tuple found in the dcc.Store
        - layout: Saved ppi graph layout in the dcc.Store

        Returns:
        - ('ppi-module-graph', 'elements'): List of elements of the ppi graph
        - ('ppi-module-graph', 'layout'): Saved ppi graph layout; otherwise 'circle' for default value
        - ('ppi-module-graph', 'style'): {'visibility': 'hidden'} for hiding the ppi graph
        - ('ppi-graph-container', 'style'): {'visibility': 'hidden'} for hiding the ppi graph container
        """

        if ppi_is_submitted:
            if submitted_algo and submitted_algo in submitted_parameter_slider:
                parameters = submitted_parameter_slider[submitted_algo]["value"]
                if not layout:
                    layout = "circle"

                return load_module_graph(
                    combined_gene_ids,
                    None,
                    submitted_network,
                    submitted_algo,
                    parameters,
                    layout,
                ) + ({"visibility": "hidden"},)

        raise PreventUpdate

    @app.callback(
        Output("ppi-download-graph-to-json", "data"),
        Input("ppi-export-graph", "n_clicks"),
        State("ppi-submitted-network", "data"),
        State("ppi-submitted-clustering-algo", "data"),
        State("ppi-submitted-parameter-slider", "data"),
        State("ppi-modules", "value"),
    )
    def download_ppi_graph_to_tsv(
        download_n_clicks,
        submitted_network,
        submitted_algo,
        submitted_parameter_slider,
        module,
    ):
        """
        Export the ppi graph in csv / tsv file format

        Parameters:
        - download_n_clicks: Number of clicks pressed on the export ppi table button
        - submitted_network: Saved ppi network found in the dcc.Store
        - submitted_algo: Saved ppi clustering algorithm found in the dcc.Store
        - submitted_parameter_slider: Saved parameter slider tuple found in the dcc.Store
        - module: Selected ppi module

        Returns:
        - ('ppi-download-graph-to-json', 'data'): PPI graph in csv / tsv file format data
        """

        if download_n_clicks >= 1:
            parameters = submitted_parameter_slider[submitted_algo]["value"]
            module_idx = int(module.split(" ")[1])
            df = pd.read_csv(
                f"{Constants.TEMP}/{submitted_network}/{submitted_algo}/modules/{parameters}/module-{module_idx}.tsv",
                sep="\t",
            )
            return dcc.send_data_frame(
                df.to_csv,
                f"[{module}] PPI Network Analysis Graph.tsv",
                index=False,
                sep="\t",
            )

        raise PreventUpdate

    # =================
    # Session-related
    # =================

    @app.callback(
        Output("ppi-graph-active-layout", "data", allow_duplicate=True),
        Output("ppi-pathway-active-tab", "data", allow_duplicate=True),
        Output("ppi-submitted-module", "data", allow_duplicate=True),
        Input("ppi-modules", "value"),
        Input("ppi-graph-layout", "value"),
        Input("ppi-modules-pathway", "active_tab"),
        State("homepage-is-submitted", "data"),
        prevent_initial_call=True,
    )
    def set_submitted_ppi_session_state(
        module, layout, active_tab, homepage_is_submitted
    ):
        """
        Sets the submitted ppi related dcc.Store variables data

        Parameters:
        - module: Selected ppi module
        - layout: Selected ppi graph layout
        - active_tab: Selected tab for ppi table
        - homepage_is_submitted: [PPI] Saved boolean value of submitted valid input

        Returns:
        - ('ppi-graph-active-layout', 'data'): Selected graph layout
        - ('ppi-pathway-active-tab', 'data'): Selected active tab for ppi table
        - ('ppi-submitted-module', 'data'): Selected ppi module
        """

        if homepage_is_submitted:
            return layout, active_tab, module

        raise PreventUpdate

    @app.callback(
        Output("ppi-graph-layout", "value"),
        Output("ppi-modules-pathway", "active_tab"),
        Input("ppi-submitted-network", "data"),
        Input("ppi-submitted-clustering-algo", "data"),
        State("ppi-is-submitted", "data"),
        State("ppi-graph-active-layout", "data"),
        State("ppi-pathway-active-tab", "data"),
    )
    def get_submitted_ppi_session_state(
        submitted_network, submitted_algo, ppi_is_submitted, layout, active_tab
    ):
        """
        Gets the [Results container] ppi related dcc.Store data and displays them

        Parameters:
        - submitted_network: Saved ppi network found in the dcc.Store
        - submitted_algo: Saved clustering algorithm found in the dcc.Store
        - layout: Saved ppi graph layout found in the dcc.Store
        - active_tab: Saved ppi active tab for the ppi table found in the dcc.Store

        Returns:
        - ('ppi-graph-layout', 'value'): Saved layout found in the dcc.Store; otherwise 'circle' for default value
        - ('ppi-modules-pathway', 'value'): Saved ppi module pathway found in the dcc.Store; otherwise 'tab-0' for default value
        """

        if ppi_is_submitted:
            if not layout:
                layout = "circle"

            if not active_tab:
                active_tab = "tab-0"

            return layout, active_tab

        raise PreventUpdate

    @app.callback(
        Output("ppi-clustering-algo", "value"),
        Output("ppi-addl-genes", "value"),
        Output("ppi-network", "value"),
        State("ppi-submitted-clustering-algo", "data"),
        State("ppi-submitted-addl-genes", "data"),
        State("ppi-submitted-network", "data"),
        Input("ppi-is-submitted", "data"),
    )
    def get_input_ppi_session_state(algo, genes, network, *_):
        """
        Gets the [Input container] ppi related dcc.Store data and displays them

        Parameters:
        - algo: Saved clustering algorithm found in the dcc.Store
        - genes: Saved ppi genes found in the dcc.Store
        - network: Saved ppi network found in the dcc.Store
        - *_: Other inputs in facilitating the saved state of the ppi input

        Returns:
        - ('ppi-clustering-algo', 'value'): Saved clustering algorithm found in the dcc.Store; otherwise 'clusterone' for default value
        - ('ppi-addl-genes', 'value'): Saved ppi additional genes found in the dcc.Store; otherwise '' for default value
        - ('ppi-network', 'value'): Saved ppi network found in the dcc.Store; otherwise 'STRING' for default value
        """

        if not algo:
            algo = "clusterone"

        if not genes:
            genes = ""

        if not network:
            network = "STRING-Physical"

        return algo, genes, network

    # =========================================
    # Network statistics / hub gene identification
    # =========================================

    @app.callback(
        Output("ppi-network-stats", "children"),
        Output("ppi-hub-genes", "children"),
        Output("ppi-network-alert", "style"),
        Output("ppi-network-alert", "children"),
        State("ppi-submitted-network", "data"),
        Input("ppi-combined-genes", "data"),
        State("ppi-is-submitted", "data"),
    )
    def display_ppi_network_stats(submitted_network, combined_gene_ids, ppi_is_submitted):
        """
        Displays topological summary statistics (node count, edge count,
        density, average node degree) and the ranked hub genes for the
        PPI subnetwork induced by the query genes.

        Also surfaces user-facing alerts, rather than an unhandled error,
        when the query gene list is empty, when none of the query genes
        are recognized in the selected PPI network, or when the query
        returns no interactions.

        Parameters:
        - submitted_network: Saved PPI network found in the dcc.Store
        - combined_gene_ids: Saved PPI combined gene ids found in the dcc.Store
        - ppi_is_submitted: [PPI] Saved boolean value of submitted valid input

        Returns:
        - ('ppi-network-stats', 'children'): Topological summary statistics
        - ('ppi-hub-genes', 'children'): Ranked list of hub genes
        - ('ppi-network-alert', 'style'): {'display': 'block'} for displaying a warning; otherwise {'display': 'none'}
        - ('ppi-network-alert', 'children'): Warning message, if any
        """

        if not ppi_is_submitted:
            raise PreventUpdate

        if not submitted_network:
            submitted_network = "STRING-Physical"

        summary = get_query_network_summary(submitted_network, combined_gene_ids or [])
        stats = summary["stats"]
        hub_genes = summary["hub_genes"]

        stats_display = [
            html.B("Nodes: "),
            f"{stats['num_nodes']:,}",
            html.Br(),
            html.B("Edges: "),
            f"{stats['num_edges']:,}",
            html.Br(),
            html.B("Network Density: "),
            f"{stats['density']:.4f}",
            html.Br(),
            html.B("Average Node Degree: "),
            f"{stats['avg_degree']:.2f}",
        ]

        if hub_genes:
            hub_genes_display = html.Ol(
                [
                    html.Li(f"{gene} (degree = {degree})")
                    for gene, degree in hub_genes
                ]
            )
        else:
            hub_genes_display = html.Span("No hub genes to display.")

        if not combined_gene_ids:
            return (
                stats_display,
                hub_genes_display,
                {"display": "block"},
                "No query genes were provided. Enter a genomic interval, or additional "
                "genes, to build a PPI network.",
            )

        if summary["is_empty"]:
            if summary["unrecognized_genes"]:
                message = (
                    "None of the submitted genes were recognized in the "
                    f"{get_user_facing_network(submitted_network)}, or they had no "
                    "recorded interactions."
                )
            else:
                message = (
                    "The submitted genes were recognized, but no interactions were "
                    f"found among them in the {get_user_facing_network(submitted_network)}."
                )

            return stats_display, hub_genes_display, {"display": "block"}, message

        if summary["unrecognized_genes"]:
            unrecognized = ", ".join(sorted(summary["unrecognized_genes"]))
            return (
                stats_display,
                hub_genes_display,
                {"display": "block"},
                f"The following gene(s) were not recognized in the "
                f"{get_user_facing_network(submitted_network)} and were excluded: "
                f"{unrecognized}.",
            )

        return stats_display, hub_genes_display, {"display": "none"}, None

    # =========================================
    # Cross-module gene handoff (from Summary)
    # =========================================

    @app.callback(
        Output("ppi-addl-genes", "value", allow_duplicate=True),
        State("ppi-carried-over-genes", "data"),
        Input("homepage-is-submitted", "data"),
        prevent_initial_call=True,
    )
    def prepopulate_ppi_addl_genes_from_summary(carried_over_genes, homepage_is_submitted):
        """
        Detects a gene list prioritized upstream on the Summary page (stored
        client-side in the 'ppi-carried-over-genes' dcc.Store) and, if present,
        pre-populates the PPI module's additional-genes input field on page load.

        Parameters:
        - carried_over_genes: Gene list stored by the Summary page's
          "Analyze in PPI Network" action
        - homepage_is_submitted: [Homepage] Saved boolean value of submitted valid input

        Returns:
        - ('ppi-addl-genes', 'value'): Semicolon-separated gene list to prepopulate
          the additional-genes textarea with
        """

        if homepage_is_submitted and carried_over_genes:
            return "; ".join(carried_over_genes)

        raise PreventUpdate

    @app.callback(
        Output("ppi-carried-over-genes", "data"),
        Output("current-analysis-page-nav", "data", allow_duplicate=True),
        Input("summary-send-genes-to-ppi", "n_clicks"),
        State("summary-results-table", "data"),
        prevent_initial_call=True,
    )
    def send_summary_genes_to_ppi(n_clicks, summary_table_data):
        """
        Handles the "Analyze in PPI Network" action on the Summary page: takes
        the gene list currently prioritized/displayed on the Summary page,
        stores it client-side, and switches the active analysis page to PPI
        Network Analysis, where it is picked up and used to pre-populate the
        additional-genes input (see prepopulate_ppi_addl_genes_from_summary).

        Parameters:
        - n_clicks: Number of clicks pressed on the "Analyze in PPI Network" button
        - summary_table_data: Rows currently displayed in the Summary results table

        Returns:
        - ('ppi-carried-over-genes', 'data'): Prioritized gene list
        - ('current-analysis-page-nav', 'data'): Constants.LABEL_PPI, to switch to the PPI page
        """

        if n_clicks and n_clicks >= 1 and summary_table_data:
            import re

            genes = []
            for row in summary_table_data:
                cell = row.get("Gene")
                if not cell or cell == NULL_PLACEHOLDER:
                    continue

                # The 'Gene' column is rendered as a Markdown link
                # (e.g., "[LOC_Os01g01010](...)"); extract the raw MSU accession
                match = re.search(r"LOC_Os\d{2}g\d{5}", cell)
                if match:
                    genes.append(match.group(0))

            return genes, Constants.LABEL_PPI

        raise PreventUpdate