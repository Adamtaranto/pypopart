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
