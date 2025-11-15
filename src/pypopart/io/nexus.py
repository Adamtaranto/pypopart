"""
NEXUS file format reader and writer for PyPopART.

Supports PopART-style NEXUS files with traits blocks.
"""

import gzip
import re
from pathlib import Path
from typing import Optional, Dict, List, Union, TextIO, Tuple

from pypopart.core.sequence import Sequence
from pypopart.core.alignment import Alignment


class NexusReader:
    """
    Reader for NEXUS format files.
    
    Supports PopART-style NEXUS with traits blocks.
    """
    
    def __init__(self, filepath: Union[str, Path], validate: bool = True):
        """
        Initialize NEXUS reader.
        
        Args:
            filepath: Path to NEXUS file
            validate: Whether to validate sequences and alignment
        """
        self.filepath = Path(filepath)
        self.validate = validate
        
        if not self.filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        self.traits = {}
    
    def _open_file(self) -> TextIO:
        """Open file handling gzip compression."""
        if self.filepath.suffix == '.gz':
            return gzip.open(self.filepath, 'rt')
        else:
            return open(self.filepath, 'r')
    
    def _parse_dimensions(self, content: str) -> Tuple[int, int]:
        """
        Parse DIMENSIONS block.
        
        Args:
            content: NEXUS file content
            
        Returns:
            Tuple of (ntax, nchar)
        """
        dimensions_match = re.search(
            r'DIMENSIONS\s+NTAX=(\d+)\s+NCHAR=(\d+)',
            content,
            re.IGNORECASE
        )
        
        if dimensions_match:
            ntax = int(dimensions_match.group(1))
            nchar = int(dimensions_match.group(2))
            return ntax, nchar
        
        return 0, 0
    
    def _parse_matrix(self, content: str) -> Dict[str, str]:
        """
        Parse MATRIX block.
        
        Args:
            content: NEXUS file content
            
        Returns:
            Dictionary mapping sequence IDs to sequence data
        """
        sequences = {}
        
        # Find MATRIX block
        matrix_match = re.search(
            r'MATRIX\s+(.*?);\s*END;',
            content,
            re.IGNORECASE | re.DOTALL
        )
        
        if not matrix_match:
            raise ValueError("No MATRIX block found in NEXUS file")
        
        matrix_content = matrix_match.group(1)
        
        # Parse sequences (handle both interleaved and sequential formats)
        current_id = None
        current_seq = []
        
        for line in matrix_content.split('\n'):
            line = line.strip()
            if not line:
                continue
            
            # Check if line starts with sequence ID
            parts = line.split(None, 1)
            if len(parts) == 2:
                seq_id, seq_data = parts
                seq_data = seq_data.replace(' ', '').replace('\t', '')
                
                if seq_id in sequences:
                    # Interleaved format - append to existing sequence
                    sequences[seq_id] += seq_data
                else:
                    # New sequence
                    sequences[seq_id] = seq_data
            elif len(parts) == 1:
                # Continuation of previous sequence
                seq_data = parts[0].replace(' ', '').replace('\t', '')
                if current_id:
                    sequences[current_id] += seq_data
        
        return sequences
    
    def _parse_traits(self, content: str) -> Dict[str, Dict[str, str]]:
        """
        Parse TRAITS block (PopART extension).
        
        Args:
            content: NEXUS file content
            
        Returns:
            Dictionary mapping sequence IDs to trait dictionaries
        """
        traits = {}
        
        # Find TRAITS block
        traits_match = re.search(
            r'BEGIN TRAITS;(.*?)END;',
            content,
            re.IGNORECASE | re.DOTALL
        )
        
        if not traits_match:
            return traits
        
        traits_content = traits_match.group(1)
        
        # Parse TRAITLABELS
        labels_match = re.search(
            r'TRAITLABELS\s+(.*?);',
            traits_content,
            re.IGNORECASE
        )
        
        trait_labels = []
        if labels_match:
            trait_labels = labels_match.group(1).split()
        
        # Parse MATRIX
        matrix_match = re.search(
            r'MATRIX\s+(.*?);',
            traits_content,
            re.IGNORECASE | re.DOTALL
        )
        
        if matrix_match:
            matrix_content = matrix_match.group(1)
            for line in matrix_content.split('\n'):
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split()
                if len(parts) >= 2:
                    seq_id = parts[0]
                    trait_values = parts[1:]
                    
                    traits[seq_id] = {}
                    for i, value in enumerate(trait_values):
                        label = trait_labels[i] if i < len(trait_labels) else f"trait_{i}"
                        traits[seq_id][label] = value
        
        return traits
    
    def read_alignment(self, progress_callback=None) -> Alignment:
        """
        Read alignment from NEXUS file.
        
        Args:
            progress_callback: Optional callback function(current, total)
            
        Returns:
            Alignment object with metadata
        """
        with self._open_file() as handle:
            content = handle.read()
        
        # Parse dimensions
        ntax, nchar = self._parse_dimensions(content)
        
        # Parse sequences
        seq_dict = self._parse_matrix(content)
        
        # Parse traits
        self.traits = self._parse_traits(content)
        
        # Create Sequence objects
        sequences = []
        count = 0
        
        for seq_id, seq_data in seq_dict.items():
            metadata = self.traits.get(seq_id, {})
            
            seq = Sequence(
                id=seq_id,
                data=seq_data.upper(),
                metadata=metadata
            )
            
            # Validation happens automatically in Sequence.__init__
            
            sequences.append(seq)
            
            count += 1
            if progress_callback:
                progress_callback(count, ntax)
        
        alignment = Alignment(sequences)
        
        if self.validate:
            alignment.validate()
        
        return alignment
    
    def get_traits(self) -> Dict[str, Dict[str, str]]:
        """
        Get parsed traits/metadata.
        
        Returns:
            Dictionary mapping sequence IDs to trait dictionaries
        """
        return self.traits


class NexusWriter:
    """
    Writer for NEXUS format files.
    
    Supports PopART-style NEXUS with traits blocks.
    """
    
    def __init__(
        self,
        filepath: Union[str, Path],
        interleaved: bool = False,
        compress: Optional[str] = None
    ):
        """
        Initialize NEXUS writer.
        
        Args:
            filepath: Output file path
            interleaved: Whether to write in interleaved format
            compress: Compression format ('gzip' or None)
        """
        self.filepath = Path(filepath)
        self.interleaved = interleaved
        self.compress = compress
        
        if compress == 'gzip' and not str(self.filepath).endswith('.gz'):
            self.filepath = Path(str(self.filepath) + '.gz')
    
    def _open_file(self) -> TextIO:
        """Open file for writing with optional compression."""
        if self.compress == 'gzip':
            return gzip.open(self.filepath, 'wt')
        else:
            return open(self.filepath, 'w')
    
    def write_alignment(
        self,
        alignment: Alignment,
        include_traits: bool = True,
        progress_callback=None
    ) -> None:
        """
        Write alignment to NEXUS file.
        
        Args:
            alignment: Alignment object
            include_traits: Whether to include traits block
            progress_callback: Optional callback function(current, total)
        """
        with self._open_file() as handle:
            # Write header
            handle.write("#NEXUS\n\n")
            
            # Write DATA block
            handle.write("BEGIN DATA;\n")
            handle.write(f"  DIMENSIONS NTAX={len(alignment)} NCHAR={alignment.length};\n")
            handle.write("  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
            handle.write("  MATRIX\n")
            
            # Write sequences
            count = 0
            for seq in alignment:
                # Pad ID to align sequences
                padded_id = seq.id.ljust(20)
                handle.write(f"    {padded_id} {seq.data}\n")
                
                count += 1
                if progress_callback:
                    progress_callback(count, len(alignment))
            
            handle.write("  ;\n")
            handle.write("END;\n")
            
            # Write TRAITS block if requested and metadata exists
            if include_traits:
                # Collect all trait keys
                all_traits = set()
                for seq in alignment:
                    all_traits.update(seq.metadata.keys())
                
                if all_traits:
                    handle.write("\n")
                    handle.write("BEGIN TRAITS;\n")
                    handle.write("  DIMENSIONS NTRAITS={};\n".format(len(all_traits)))
                    
                    # Write trait labels
                    trait_list = sorted(all_traits)
                    handle.write("  TRAITLABELS {};\n".format(' '.join(trait_list)))
                    
                    # Write trait matrix
                    handle.write("  MATRIX\n")
                    for seq in alignment:
                        padded_id = seq.id.ljust(20)
                        trait_values = [
                            str(seq.metadata.get(trait, '?'))
                            for trait in trait_list
                        ]
                        handle.write(f"    {padded_id} {' '.join(trait_values)}\n")
                    
                    handle.write("  ;\n")
                    handle.write("END;\n")
