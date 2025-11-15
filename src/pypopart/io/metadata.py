"""
Metadata file reader and writer for PyPopART.

Handles CSV-based metadata/traits files that can be linked to sequences.
"""

import csv
import gzip
from pathlib import Path
from typing import Dict, List, Optional, TextIO, Union

from pypopart.core.alignment import Alignment


class MetadataReader:
    """
    Reader for CSV-based metadata files.
    """

    def __init__(
        self,
        filepath: Union[str, Path],
        id_column: str = 'id',
        delimiter: str = ',',
        validate: bool = True,
    ):
        """
        Initialize metadata reader.

        Args:
            filepath: Path to metadata CSV file
            id_column: Name of column containing sequence IDs
            delimiter: CSV delimiter character
            validate: Whether to validate metadata
        """
        self.filepath = Path(filepath)
        self.id_column = id_column
        self.delimiter = delimiter
        self.validate = validate

        if not self.filepath.exists():
            raise FileNotFoundError(f'File not found: {filepath}')

    def _open_file(self) -> TextIO:
        """Open file handling gzip compression."""
        if self.filepath.suffix == '.gz':
            return gzip.open(self.filepath, 'rt')
        else:
            return open(self.filepath, 'r', encoding='utf-8')

    def read_metadata(self) -> Dict[str, Dict[str, str]]:
        """
        Read metadata from CSV file.

        Returns:
            Dictionary mapping sequence IDs to metadata dictionaries
        """
        metadata = {}

        with self._open_file() as handle:
            reader = csv.DictReader(handle, delimiter=self.delimiter)

            if self.id_column not in reader.fieldnames:
                raise ValueError(
                    f"ID column '{self.id_column}' not found in CSV. "
                    f'Available columns: {reader.fieldnames}'
                )

            for row in reader:
                seq_id = row[self.id_column]

                # Remove ID column from metadata
                row_metadata = {k: v for k, v in row.items() if k != self.id_column}

                if seq_id in metadata and self.validate:
                    raise ValueError(f'Duplicate sequence ID in metadata: {seq_id}')

                metadata[seq_id] = row_metadata

        return metadata

    def apply_to_alignment(self, alignment: Alignment) -> None:
        """
        Apply metadata to sequences in alignment.

        Args:
            alignment: Alignment object to update
        """
        metadata = self.read_metadata()

        missing_ids = []
        for seq in alignment:
            if seq.id in metadata:
                # Update sequence metadata
                seq.metadata.update(metadata[seq.id])
            elif self.validate:
                missing_ids.append(seq.id)

        if missing_ids and self.validate:
            raise ValueError(
                f'Metadata missing for {len(missing_ids)} sequences: '
                f'{", ".join(missing_ids[:5])}{"..." if len(missing_ids) > 5 else ""}'
            )


class MetadataWriter:
    """
    Writer for CSV-based metadata files.
    """

    def __init__(
        self,
        filepath: Union[str, Path],
        id_column: str = 'id',
        delimiter: str = ',',
        compress: Optional[str] = None,
    ):
        """
        Initialize metadata writer.

        Args:
            filepath: Output file path
            id_column: Name of column for sequence IDs
            delimiter: CSV delimiter character
            compress: Compression format ('gzip' or None)
        """
        self.filepath = Path(filepath)
        self.id_column = id_column
        self.delimiter = delimiter
        self.compress = compress

        if compress == 'gzip' and not str(self.filepath).endswith('.gz'):
            self.filepath = Path(str(self.filepath) + '.gz')

    def _open_file(self) -> TextIO:
        """Open file for writing with optional compression."""
        if self.compress == 'gzip':
            return gzip.open(self.filepath, 'wt', encoding='utf-8')
        else:
            return open(self.filepath, 'w', encoding='utf-8', newline='')

    def write_metadata(
        self,
        metadata: Dict[str, Dict[str, str]],
        trait_order: Optional[List[str]] = None,
    ) -> None:
        """
        Write metadata to CSV file.

        Args:
            metadata: Dictionary mapping sequence IDs to metadata dictionaries
            trait_order: Optional list specifying order of trait columns
        """
        if not metadata:
            raise ValueError('No metadata to write')

        # Determine all trait keys
        all_traits = set()
        for traits in metadata.values():
            all_traits.update(traits.keys())

        # Order traits
        if trait_order:
            # Use specified order, then add any remaining traits
            fieldnames = [self.id_column] + trait_order
            remaining = sorted(all_traits - set(trait_order))
            fieldnames.extend(remaining)
        else:
            fieldnames = [self.id_column] + sorted(all_traits)

        with self._open_file() as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                delimiter=self.delimiter,
                extrasaction='ignore',
            )

            writer.writeheader()

            for seq_id, traits in sorted(metadata.items()):
                row = {self.id_column: seq_id}
                row.update(traits)
                writer.writerow(row)

    def write_from_alignment(
        self, alignment: Alignment, trait_order: Optional[List[str]] = None
    ) -> None:
        """
        Extract and write metadata from alignment.

        Args:
            alignment: Alignment object
            trait_order: Optional list specifying order of trait columns
        """
        metadata = {}

        for seq in alignment:
            if seq.metadata:
                metadata[seq.id] = seq.metadata

        if not metadata:
            raise ValueError('No metadata found in alignment')

        self.write_metadata(metadata, trait_order)
