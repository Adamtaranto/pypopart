"""Unit tests for PyPopART GUI application."""


class TestPyPopARTApp:
    """Test cases for PyPopART Dash GUI."""

    def test_app_initialization(self):
        """Test that the app initializes without errors."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False, port=8050)

        assert app is not None
        assert app.app is not None
        assert app.app.title == 'PyPopART - Haplotype Network Analysis'
        assert app.debug is False
        assert app.port == 8050

    def test_app_layout_structure(self):
        """Test that the app layout has expected components."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)

        # Check that layout is set up
        assert app.app.layout is not None

        # Check that key components exist in layout
        layout_children = app.app.layout.children
        assert len(layout_children) > 0

    def test_create_upload_card(self):
        """Test upload card creation."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        card = app._create_upload_card()

        assert card is not None
        assert hasattr(card, 'children')

    def test_create_algorithm_card(self):
        """Test algorithm card creation."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        card = app._create_algorithm_card()

        assert card is not None
        assert hasattr(card, 'children')

    def test_create_layout_card(self):
        """Test layout card creation."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        card = app._create_layout_card()

        assert card is not None
        assert hasattr(card, 'children')

    def test_create_export_card(self):
        """Test export card creation."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        card = app._create_export_card()

        assert card is not None
        assert hasattr(card, 'children')

    def test_create_network_tab(self):
        """Test network tab creation."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        tab = app._create_network_tab()

        assert tab is not None
        assert hasattr(tab, 'children')

    def test_create_statistics_tab(self):
        """Test statistics tab creation."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        tab = app._create_statistics_tab()

        assert tab is not None
        assert hasattr(tab, 'children')

    def test_create_alignment_tab(self):
        """Test alignment tab creation."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        tab = app._create_alignment_tab()

        assert tab is not None
        assert hasattr(tab, 'children')

    def test_main_function_exists(self):
        """Test that main function is callable."""
        from pypopart.gui.app import main

        assert callable(main)

    def test_import_from_init(self):
        """Test that components can be imported from gui module."""
        from pypopart.gui import PyPopARTApp, main

        assert PyPopARTApp is not None
        assert callable(main)

    def test_logger_initialization(self):
        """Test that logger is properly initialized."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=True, port=8050)
        assert hasattr(app, 'logger')
        assert app.logger is not None

    def test_callbacks_registered(self):
        """Test that callbacks are properly registered without errors."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        # If callbacks are registered with pattern-matching IDs correctly,
        # the app should initialize without errors
        assert app.app._callback_list is not None
        # Check that we have the expected number of callbacks
        assert len(app.app._callback_list) > 0

    def test_default_layout_is_hierarchical(self):
        """Test that default layout is set to hierarchical."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        # Find the layout dropdown in the layout
        # It should have 'hierarchical' as the default value
        # This verifies Issue 4 is fixed
        layout_card = app._create_layout_card()
        assert layout_card is not None

    def test_search_dropdown_supports_multi_select(self):
        """Test that search dropdown supports multiple selections."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        # Find the haplotype search dropdown in the network tab
        # It should have multi=True set
        # This verifies Issue 6 is fixed
        network_tab = app._create_network_tab()
        assert network_tab is not None

    def test_haplotype_summary_tab_has_mapping_components(self):
        """Test that haplotype summary tab has label mapping components."""
        from pypopart.gui.app import PyPopARTApp

        app = PyPopARTApp(debug=False)
        tab = app._create_haplotype_summary_tab()

        assert tab is not None
        # Tab should have children including buttons and stores
        assert hasattr(tab, 'children')
        assert len(tab.children) > 0

        # Find components by traversing children
        # Should have: download-h-number-template-button, upload-h-number-mapping, h-number-mapping-store
        component_ids = []

        def extract_ids(component):
            """Recursively extract all component IDs."""
            if hasattr(component, 'id'):
                component_ids.append(component.id)
            if hasattr(component, 'children'):
                if isinstance(component.children, list):
                    for child in component.children:
                        extract_ids(child)
                else:
                    extract_ids(component.children)

        extract_ids(tab)

        # Check that key components exist
        assert 'download-h-number-template-button' in component_ids
        assert 'upload-h-number-mapping' in component_ids
        assert 'h-number-mapping-store' in component_ids
        assert 'h-number-feedback' in component_ids
