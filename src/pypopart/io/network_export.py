"""
Network export module for PyPopART.

Provides exporters for various network file formats.
"""

import json
import csv
from pathlib import Path
from typing import Optional, Union, Dict, Any
import networkx as nx

from pypopart.core.graph import HaplotypeNetwork


class GraphMLExporter:
    """
    Export haplotype networks to GraphML format.
    """
    
    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize GraphML exporter.
        
        Args:
            filepath: Output file path
        """
        self.filepath = Path(filepath)
    
    def export(self, network: HaplotypeNetwork) -> None:
        """
        Export network to GraphML format.
        
        Args:
            network: HaplotypeNetwork object
        """
        # Convert network to NetworkX graph if needed
        graph = network.graph if hasattr(network, 'graph') else network
        
        # Write to GraphML
        nx.write_graphml(graph, self.filepath)


class GMLExporter:
    """
    Export haplotype networks to GML format.
    """
    
    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize GML exporter.
        
        Args:
            filepath: Output file path
        """
        self.filepath = Path(filepath)
    
    def export(self, network: HaplotypeNetwork) -> None:
        """
        Export network to GML format.
        
        Args:
            network: HaplotypeNetwork object
        """
        graph = network.graph if hasattr(network, 'graph') else network
        
        # Write to GML
        nx.write_gml(graph, self.filepath)


class CytoscapeExporter:
    """
    Export haplotype networks to Cytoscape JSON format.
    """
    
    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize Cytoscape exporter.
        
        Args:
            filepath: Output file path
        """
        self.filepath = Path(filepath)
    
    def export(self, network: HaplotypeNetwork) -> None:
        """
        Export network to Cytoscape JSON format.
        
        Args:
            network: HaplotypeNetwork object
        """
        graph = network.graph if hasattr(network, 'graph') else network
        
        # Convert to Cytoscape JSON format
        cytoscape_data = nx.cytoscape_data(graph)
        
        with open(self.filepath, 'w') as f:
            json.dump(cytoscape_data, f, indent=2)


class JSONExporter:
    """
    Export haplotype networks to JSON format.
    """
    
    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize JSON exporter.
        
        Args:
            filepath: Output file path
        """
        self.filepath = Path(filepath)
    
    def export(
        self,
        network: HaplotypeNetwork,
        include_layout: bool = True
    ) -> None:
        """
        Export network to JSON format.
        
        Args:
            network: HaplotypeNetwork object
            include_layout: Whether to include node layout positions
        """
        graph = network.graph if hasattr(network, 'graph') else network
        
        # Build JSON structure
        data = {
            'nodes': [],
            'edges': [],
            'metadata': {}
        }
        
        # Add nodes
        for node, attrs in graph.nodes(data=True):
            node_data = {
                'id': str(node),
                'attributes': dict(attrs)
            }
            data['nodes'].append(node_data)
        
        # Add edges
        for source, target, attrs in graph.edges(data=True):
            edge_data = {
                'source': str(source),
                'target': str(target),
                'attributes': dict(attrs)
            }
            data['edges'].append(edge_data)
        
        # Add graph metadata
        if hasattr(graph, 'graph'):
            data['metadata'] = dict(graph.graph)
        
        with open(self.filepath, 'w') as f:
            json.dump(data, f, indent=2)


class CSVExporter:
    """
    Export haplotype network statistics to CSV format.
    """
    
    def __init__(self, filepath: Union[str, Path]):
        """
        Initialize CSV exporter.
        
        Args:
            filepath: Output file path
        """
        self.filepath = Path(filepath)
    
    def export_nodes(self, network: HaplotypeNetwork) -> None:
        """
        Export node attributes to CSV.
        
        Args:
            network: HaplotypeNetwork object
        """
        graph = network.graph if hasattr(network, 'graph') else network
        
        # Collect all attribute keys
        all_keys = set()
        for node, attrs in graph.nodes(data=True):
            all_keys.update(attrs.keys())
        
        fieldnames = ['node_id'] + sorted(all_keys)
        
        with open(self.filepath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for node, attrs in graph.nodes(data=True):
                row = {'node_id': str(node)}
                row.update(attrs)
                writer.writerow(row)
    
    def export_edges(self, network: HaplotypeNetwork) -> None:
        """
        Export edge attributes to CSV.
        
        Args:
            network: HaplotypeNetwork object
        """
        graph = network.graph if hasattr(network, 'graph') else network
        
        # Collect all attribute keys
        all_keys = set()
        for source, target, attrs in graph.edges(data=True):
            all_keys.update(attrs.keys())
        
        fieldnames = ['source', 'target'] + sorted(all_keys)
        
        with open(self.filepath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for source, target, attrs in graph.edges(data=True):
                row = {'source': str(source), 'target': str(target)}
                row.update(attrs)
                writer.writerow(row)
    
    def export_statistics(
        self,
        network: HaplotypeNetwork,
        statistics: Dict[str, Any]
    ) -> None:
        """
        Export network statistics to CSV.
        
        Args:
            network: HaplotypeNetwork object
            statistics: Dictionary of statistics to export
        """
        with open(self.filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Statistic', 'Value'])
            
            for key, value in sorted(statistics.items()):
                writer.writerow([key, value])
