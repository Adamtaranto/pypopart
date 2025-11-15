"""
Unit tests for geographic layout algorithms.
"""

import pytest

from pypopart.core.graph import HaplotypeNetwork
from pypopart.core.haplotype import Haplotype
from pypopart.core.sequence import Sequence
from pypopart.layout.algorithms import GeographicLayout


class TestGeographicLayout:
    """Test geographic layout algorithm."""
    
    @pytest.fixture
    def simple_network(self):
        """Create a simple test network."""
        network = HaplotypeNetwork()
        
        # Create haplotypes
        seq_a = Sequence('A', 'ATCGATCG')
        hap_a = Haplotype(seq_a, sample_ids=['s1', 's2', 's3', 's4', 's5'])
        
        seq_b = Sequence('B', 'ATCGATCG')
        hap_b = Haplotype(seq_b, sample_ids=['s6', 's7', 's8'])
        
        seq_c = Sequence('C', 'ATCGATCG')
        hap_c = Haplotype(seq_c, sample_ids=['s9', 's10'])
        
        # Add to network
        network.add_haplotype(hap_a)
        network.add_haplotype(hap_b)
        network.add_haplotype(hap_c)
        
        # Add edges
        network.add_edge('A', 'B', mutations=1)
        network.add_edge('B', 'C', mutations=2)
        
        return network
    
    def test_compute_with_coordinates(self, simple_network):
        """Test layout computation with provided coordinates."""
        coordinates = {
            'A': (40.7128, -74.0060),  # New York
            'B': (34.0522, -118.2437),  # Los Angeles
            'C': (51.5074, -0.1278),  # London
        }
        
        geo_layout = GeographicLayout(simple_network)
        layout = geo_layout.compute(coordinates=coordinates, projection='platecarree')
        
        assert len(layout) == 3
        assert 'A' in layout
        assert 'B' in layout
        assert 'C' in layout
        
        # Verify positions are tuples
        for pos in layout.values():
            assert isinstance(pos, tuple)
            assert len(pos) == 2
    
    def test_compute_mercator_projection(self, simple_network):
        """Test Mercator projection."""
        coordinates = {
            'A': (40.7128, -74.0060),
            'B': (34.0522, -118.2437),
        }
        
        geo_layout = GeographicLayout(simple_network)
        layout = geo_layout.compute(coordinates=coordinates, projection='mercator')
        
        # Mercator should transform y coordinates
        assert layout['A'][0] == -74.0060  # x (longitude) unchanged
        assert layout['A'][1] != 40.7128  # y (latitude) transformed
    
    def test_compute_platecarree_projection(self, simple_network):
        """Test PlateCarree (equirectangular) projection."""
        coordinates = {
            'A': (40.7128, -74.0060),
            'B': (34.0522, -118.2437),
        }
        
        geo_layout = GeographicLayout(simple_network)
        layout = geo_layout.compute(coordinates=coordinates, projection='platecarree')
        
        # PlateCarree should not transform coordinates
        assert layout['A'] == (-74.0060, 40.7128)
        assert layout['B'] == (-118.2437, 34.0522)
    
    def test_compute_missing_coordinates(self, simple_network):
        """Test layout with missing coordinates for some nodes."""
        coordinates = {
            'A': (40.7128, -74.0060),
            'B': (34.0522, -118.2437),
            # C has no coordinates
        }
        
        geo_layout = GeographicLayout(simple_network)
        layout = geo_layout.compute(coordinates=coordinates, projection='platecarree')
        
        # All nodes should have positions (C gets a default position)
        assert len(layout) == 3
        assert 'A' in layout
        assert 'B' in layout
        assert 'C' in layout
    
    def test_compute_no_coordinates_error(self, simple_network):
        """Test error when no coordinates are provided."""
        geo_layout = GeographicLayout(simple_network)
        
        with pytest.raises(ValueError, match='No geographic coordinates'):
            geo_layout.compute()
    
    def test_invalid_latitude(self, simple_network):
        """Test error handling for invalid latitude."""
        coordinates = {
            'A': (95.0, -74.0060),  # Invalid latitude
        }
        
        geo_layout = GeographicLayout(simple_network)
        
        with pytest.raises(ValueError, match='Invalid latitude'):
            geo_layout.compute(coordinates=coordinates)
    
    def test_invalid_longitude(self, simple_network):
        """Test error handling for invalid longitude."""
        coordinates = {
            'A': (40.7128, -200.0),  # Invalid longitude
        }
        
        geo_layout = GeographicLayout(simple_network)
        
        with pytest.raises(ValueError, match='Invalid longitude'):
            geo_layout.compute(coordinates=coordinates)
    
    def test_extract_coordinates_from_metadata(self, simple_network):
        """Test extraction of coordinates from node metadata."""
        # Add coordinates to node metadata as strings (as they would be from CSV)
        simple_network._graph.nodes['A']['latitude'] = '40.7128'
        simple_network._graph.nodes['A']['longitude'] = '-74.0060'
        simple_network._graph.nodes['B']['latitude'] = '34.0522'
        simple_network._graph.nodes['B']['longitude'] = '-118.2437'
        
        geo_layout = GeographicLayout(simple_network)
        coordinates = geo_layout._extract_coordinates_from_metadata()
        
        assert len(coordinates) == 2
        assert coordinates['A'] == (40.7128, -74.0060)
        assert coordinates['B'] == (34.0522, -118.2437)
    
    def test_jitter_amount(self, simple_network):
        """Test jitter functionality."""
        coordinates = {
            'A': (40.0, -74.0),
            'B': (40.0, -74.0),  # Same location as A
        }
        
        geo_layout = GeographicLayout(simple_network)
        
        # Without jitter, positions should be exactly the same
        layout1 = geo_layout.compute(
            coordinates=coordinates, 
            projection='platecarree',
            jitter_amount=0.0
        )
        assert layout1['A'] == layout1['B']
        
        # With jitter, positions should differ
        # Note: This is probabilistic, but very unlikely to be exactly the same
        layout2 = geo_layout.compute(
            coordinates=coordinates, 
            projection='platecarree',
            jitter_amount=0.5
        )
        # Positions should be close but not identical
        assert abs(layout2['A'][0] - layout2['B'][0]) < 1.0
        assert abs(layout2['A'][1] - layout2['B'][1]) < 1.0
