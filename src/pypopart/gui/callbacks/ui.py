"""Chrome-level UI behaviour: collapsing and restoring the control panel."""

from dash import Input, Output, State


def register(app, logger) -> None:
    """
    Register UI chrome callbacks on the Dash app.

    Parameters
    ----------
    app : dash.Dash
        The Dash application.
    logger : logging.Logger
        Application logger.
    """
    # Clientside on purpose. The sidebar carries CSS 'resize: horizontal',
    # so a drag writes an *inline* width straight onto the DOM node. A
    # server callback returning a new style dict would make React re-apply
    # the prop and wipe that dragged width on every toggle. Touching
    # el.style here leaves Dash's prop untouched, so the width survives.
    app.clientside_callback(
        """
        function(nClicks, collapsed) {
            if (!nClicks) {
                return window.dash_clientside.no_update;
            }

            const el = document.getElementById('sidebar-panel');
            if (!el) {
                return window.dash_clientside.no_update;
            }

            const nowCollapsed = !collapsed;
            if (nowCollapsed) {
                // Stash whatever width the user dragged to, then hide.
                window.pypopartSidebarWidth =
                    el.style.width || el.getBoundingClientRect().width + 'px';
                el.style.display = 'none';
            } else {
                el.style.display = '';
                if (window.pypopartSidebarWidth) {
                    el.style.width = window.pypopartSidebarWidth;
                }
            }

            // The debounced window-resize handler installed alongside the
            // network already calls cy.resize()/cy.fit(), so reuse it.
            setTimeout(function() {
                window.dispatchEvent(new Event('resize'));
            }, 0);

            // That handler only exists once a network has been drawn, and
            // it is debounced, so also nudge the instance directly.
            try {
                const cy = document.getElementById('network-graph')._cyreg.cy;
                if (cy) {
                    cy.resize();
                    cy.fit(null, 50);
                }
            } catch (e) {
                // No network yet; nothing to re-fit.
            }

            return [
                nowCollapsed,
                nowCollapsed ? '\\u00bb' : '\\u00ab',
                nowCollapsed ? 'Show the control panel' : 'Hide the control panel',
            ];
        }
        """,
        [
            Output('sidebar-collapsed', 'data'),
            Output('sidebar-toggle', 'children'),
            Output('sidebar-toggle', 'title'),
        ],
        Input('sidebar-toggle', 'n_clicks'),
        State('sidebar-collapsed', 'data'),
        prevent_initial_call=True,
    )

    # Publishes the grid settings onto `window` so the Cytoscape event
    # handlers below can read them without a round trip on every drag.
    app.clientside_callback(
        """
        function(enabled, size) {
            window.pypopartGrid = {
                enabled: !!enabled,
                size: size || 0,
            };
            if (window.pypopartDrawGrid) {
                window.pypopartDrawGrid();
            }
            return window.dash_clientside.no_update;
        }
        """,
        Output('grid-size', 'className'),
        [
            Input('snap-to-grid-toggle', 'value'),
            Input('grid-size', 'value'),
        ],
    )

    # Snapping on drop, plus the grid drawn behind the network. The grid is
    # a CSS background on the Cytoscape container, so it has to be redrawn
    # on pan and zoom or it would only line up with the nodes at zoom 1.
    app.clientside_callback(
        """
        function() {
            if (window.pypopartGridSetup) {
                return window.dash_clientside.no_update;
            }
            window.pypopartGridSetup = true;

            setTimeout(function() {
                try {
                    const cy = document.getElementById('network-graph')._cyreg.cy;
                    if (!cy) { return; }
                    const container = cy.container();

                    function settings() {
                        return window.pypopartGrid || {enabled: false, size: 0};
                    }

                    window.pypopartDrawGrid = function() {
                        const grid = settings();
                        if (!grid.enabled || !grid.size) {
                            container.style.backgroundImage = '';
                            return;
                        }
                        const step = grid.size * cy.zoom();
                        const pan = cy.pan();
                        const line = 'rgba(0, 48, 73, 0.12)';
                        container.style.backgroundImage =
                            'linear-gradient(to right, ' + line + ' 1px, transparent 1px),' +
                            'linear-gradient(to bottom, ' + line + ' 1px, transparent 1px)';
                        container.style.backgroundSize = step + 'px ' + step + 'px';
                        container.style.backgroundPosition =
                            pan.x + 'px ' + pan.y + 'px';
                    };

                    cy.on('pan zoom resize', window.pypopartDrawGrid);

                    cy.on('dragfree', 'node', function(evt) {
                        const grid = settings();
                        if (!grid.enabled || !grid.size) { return; }
                        const pos = evt.target.position();
                        evt.target.position({
                            x: Math.round(pos.x / grid.size) * grid.size,
                            y: Math.round(pos.y / grid.size) * grid.size,
                        });
                    });

                    window.pypopartDrawGrid();
                } catch (e) {
                    console.log('Error setting up grid:', e);
                }
            }, 500);

            return window.dash_clientside.no_update;
        }
        """,
        Output('snap-to-grid-toggle', 'className'),
        Input('network-graph', 'elements'),
        prevent_initial_call=True,
    )
