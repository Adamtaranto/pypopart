"""
Unit tests for NEXUS file I/O.
"""

import pytest
from pathlib import Path
from pypopart.io.nexus import NexusReader, NexusWriter
from pypopart.core.sequence import Sequence
from pypopart.core.alignment import Alignment


class TestNexusReader:
    """Test cases for NexusReader."""
    
    def test_read_simple_nexus(self, tmp_path):
        """Test reading a simple NEXUS file."""
        nexus_file = tmp_path / "test.nex"
        nexus_file.write_text(
            "#NEXUS\n\n"
            "BEGIN DATA;\n"
            "  DIMENSIONS NTAX=2 NCHAR=4;\n"
            "  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n"
            "  MATRIX\n"
            "    seq1 ATCG\n"
            "    seq2 ATCC\n"
            "  ;\n"
            "END;\n"
        )
        
        reader = NexusReader(nexus_file)
        alignment = reader.read_alignment()
        
        assert len(alignment) == 2
        assert alignment.length == 4
        assert alignment[0].id == "seq1"
        assert alignment[0].data == "ATCG"
    
    def test_read_nexus_with_traits(self, tmp_path):
        """Test reading NEXUS with traits block."""
        nexus_file = tmp_path / "test.nex"
        nexus_file.write_text(
            "#NEXUS\n\n"
            "BEGIN DATA;\n"
            "  DIMENSIONS NTAX=2 NCHAR=4;\n"
            "  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n"
            "  MATRIX\n"
            "    seq1 ATCG\n"
            "    seq2 ATCC\n"
            "  ;\n"
            "END;\n\n"
            "BEGIN TRAITS;\n"
            "  DIMENSIONS NTRAITS=1;\n"
            "  TRAITLABELS population;\n"
            "  MATRIX\n"
            "    seq1 A\n"
            "    seq2 B\n"
            "  ;\n"
            "END;\n"
        )
        
        reader = NexusReader(nexus_file)
        alignment = reader.read_alignment()
        
        assert alignment[0].metadata['population'] == 'A'
        assert alignment[1].metadata['population'] == 'B'
    
    def test_read_interleaved_format(self, tmp_path):
        """Test reading interleaved NEXUS format."""
        nexus_file = tmp_path / "test.nex"
        nexus_file.write_text(
            "#NEXUS\n\n"
            "BEGIN DATA;\n"
            "  DIMENSIONS NTAX=2 NCHAR=8;\n"
            "  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n"
            "  MATRIX\n"
            "    seq1 ATCG\n"
            "    seq2 GTCC\n"
            "    seq1 ATAT\n"
            "    seq2 GCGC\n"
            "  ;\n"
            "END;\n"
        )
        
        reader = NexusReader(nexus_file)
        alignment = reader.read_alignment()
        
        assert alignment.length == 8
        assert alignment[0].data == "ATCGATAT"
        assert alignment[1].data == "GTCCGCGC"
    
    def test_file_not_found(self):
        """Test error handling for missing file."""
        with pytest.raises(FileNotFoundError):
            NexusReader("nonexistent.nex")


class TestNexusWriter:
    """Test cases for NexusWriter."""
    
    def test_write_alignment(self, tmp_path):
        """Test writing alignment to NEXUS."""
        alignment = Alignment([
            Sequence("seq1", "ATCG"),
            Sequence("seq2", "ATCC")
        ])
        
        output_file = tmp_path / "output.nex"
        writer = NexusWriter(output_file)
        writer.write_alignment(alignment)
        
        assert output_file.exists()
        
        # Read back and verify
        reader = NexusReader(output_file)
        read_alignment = reader.read_alignment()
        assert len(read_alignment) == 2
        assert read_alignment.length == 4
    
    def test_write_with_traits(self, tmp_path):
        """Test writing alignment with traits."""
        alignment = Alignment([
            Sequence("seq1", "ATCG", metadata={"pop": "A"}),
            Sequence("seq2", "ATCC", metadata={"pop": "B"})
        ])
        
        output_file = tmp_path / "output.nex"
        writer = NexusWriter(output_file)
        writer.write_alignment(alignment, include_traits=True)
        
        # Read back and verify traits
        reader = NexusReader(output_file)
        read_alignment = reader.read_alignment()
        assert read_alignment[0].metadata['pop'] == 'A'
        assert read_alignment[1].metadata['pop'] == 'B'
    
    def test_write_without_traits(self, tmp_path):
        """Test writing alignment without traits block."""
        alignment = Alignment([
            Sequence("seq1", "ATCG", metadata={"pop": "A"}),
            Sequence("seq2", "ATCC", metadata={"pop": "B"})
        ])
        
        output_file = tmp_path / "output.nex"
        writer = NexusWriter(output_file)
        writer.write_alignment(alignment, include_traits=False)
        
        content = output_file.read_text()
        assert "TRAITS" not in content
