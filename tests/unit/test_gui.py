"""
Unit tests for PyPopART GUI application.
"""



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
