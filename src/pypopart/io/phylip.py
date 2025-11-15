"""
PHYLIP file format reader and writer for PyPopART.
"""

import gzip
from pathlib import Path
from typing import Optional, TextIO, Union

from pypopart.core.alignment import Alignment
from pypopart.core.sequence import Sequence


class PhylipReader:
    """
    Reader for PHYLIP format sequence files.

    Supports both sequential and interleaved formats.
    """

    def __init__(
        self, filepath: Union[str, Path], strict: bool = False, validate: bool = True
    ):
        """
        Initialize PHYLIP reader.

        Parameters
        ----------
        filepath :
            Path to PHYLIP file.
        strict :
            Whether to use strict format (10-char IDs).
        validate :
            Whether to validate sequences and alignment.
        """
        self.filepath = Path(filepath)
        self.strict = strict
        self.validate = validate

        if not self.filepath.exists():
            raise FileNotFoundError(f'File not found: {filepath}')

    def _open_file(self) -> TextIO:
        """Open file handling gzip compression."""
        if self.filepath.suffix == '.gz':
            return gzip.open(self.filepath, 'rt')
        else:
            return open(self.filepath, 'r')

    def read_alignment(self, progress_callback=None) -> Alignment:
        """
            Read alignment from PHYLIP file.

            Parameters
            ----------
            progress_callback :
                Optional callback function(current, total).

        Returns
        -------
            Alignment object.
        """
        with self._open_file() as handle:
            lines = [line.rstrip() for line in handle if line.strip()]

        if not lines:
            raise ValueError('Empty PHYLIP file')

        # Parse header
        header_parts = lines[0].split()
        if len(header_parts) != 2:
            raise ValueError('Invalid PHYLIP header format')

        ntax = int(header_parts[0])
        nchar = int(header_parts[1])

        # Parse sequences
        sequences = {}
        current_line = 1

        # First pass: read sequence IDs and initial data
        for _i in range(ntax):
            if current_line >= len(lines):
                raise ValueError(f'Unexpected end of file (expected {ntax} sequences)')

            line = lines[current_line]

            if self.strict:
                # Strict format: first 10 characters are ID
                seq_id = line[:10].strip()
                seq_data = line[10:].replace(' ', '').replace('\t', '')
            else:
                # Relaxed format: ID separated by whitespace
                parts = line.split(None, 1)
                if len(parts) < 1:
                    raise ValueError(f'Invalid sequence line: {line}')

                seq_id = parts[0]
                seq_data = (
                    parts[1].replace(' ', '').replace('\t', '')
                    if len(parts) > 1
                    else ''
                )

            sequences[seq_id] = seq_data
            current_line += 1

        # Check if interleaved (more lines after initial block)
        if current_line < len(lines):
            # Interleaved format: read additional blocks
            seq_ids = list(sequences.keys())

            while current_line < len(lines):
                for seq_id in seq_ids:
                    if current_line >= len(lines):
                        break

                    line = lines[current_line].replace(' ', '').replace('\t', '')
                    sequences[seq_id] += line
                    current_line += 1

        # Create Sequence objects
        seq_list = []
        count = 0

        for seq_id, seq_data in sequences.items():
            if len(seq_data) != nchar:
                raise ValueError(
                    f"Sequence {seq_id} length {len(seq_data)} doesn't match expected {nchar}"
                )

            seq = Sequence(id=seq_id, data=seq_data.upper())

            # Validation happens automatically in Sequence.__init__

            seq_list.append(seq)

            count += 1
            if progress_callback:
                progress_callback(count, ntax)

        alignment = Alignment(seq_list)

        if self.validate:
            alignment.validate()

        return alignment


class PhylipWriter:
    """
    Writer for PHYLIP format sequence files.
    """

    def __init__(
        self,
        filepath: Union[str, Path],
        strict: bool = False,
        interleaved: bool = False,
        line_length: int = 60,
        compress: Optional[str] = None,
    ):
        """
        Initialize PHYLIP writer.

        Parameters
        ----------
        filepath :
            Output file path.
        strict :
            Whether to use strict format (10-char IDs).
        interleaved :
            Whether to write in interleaved format.
        line_length :
            Line length for interleaved format.
        compress :
            Compression format ('gzip' or None).
        """
        self.filepath = Path(filepath)
        self.strict = strict
        self.interleaved = interleaved
        self.line_length = line_length
        self.compress = compress

        if compress == 'gzip' and not str(self.filepath).endswith('.gz'):
            self.filepath = Path(str(self.filepath) + '.gz')

    def _open_file(self) -> TextIO:
        """Open file for writing with optional compression."""
        if self.compress == 'gzip':
            return gzip.open(self.filepath, 'wt')
        else:
            return open(self.filepath, 'w')

    def write_alignment(self, alignment: Alignment, progress_callback=None) -> None:
        """
        Write alignment to PHYLIP file.

        Parameters
        ----------
        alignment :
            Alignment object.
        progress_callback :
            Optional callback function(current, total).
        """
        with self._open_file() as handle:
            # Write header
            handle.write(f' {len(alignment)} {alignment.length}\n')

            if self.interleaved:
                # Interleaved format
                num_blocks = (
                    alignment.length + self.line_length - 1
                ) // self.line_length

                for block in range(num_blocks):
                    start = block * self.line_length
                    end = min(start + self.line_length, alignment.length)

                    for i, seq in enumerate(alignment):
                        if block == 0:
                            # First block: include ID
                            if self.strict:
                                seq_id = seq.id[:10].ljust(10)
                            else:
                                seq_id = seq.id.ljust(20)
                            handle.write(f'{seq_id} {seq.data[start:end]}\n')
                        else:
                            # Subsequent blocks: just sequence data
                            handle.write(f'{seq.data[start:end]}\n')

                        if progress_callback:
                            progress_callback(i + 1, len(alignment))

                    # Blank line between blocks
                    if block < num_blocks - 1:
                        handle.write('\n')
            else:
                # Sequential format
                for i, seq in enumerate(alignment):
                    if self.strict:
                        seq_id = seq.id[:10].ljust(10)
                    else:
                        seq_id = seq.id.ljust(20)

                    handle.write(f'{seq_id} {seq.data}\n')

                    if progress_callback:
                        progress_callback(i + 1, len(alignment))
