# Contributing to PyPopART

Thank you for your interest in contributing to PyPopART! This document provides guidelines and instructions for contributing.

## Ways to Contribute

- **Report bugs**: File issues on GitHub
- **Suggest features**: Propose new functionality
- **Fix bugs**: Submit pull requests
- **Improve documentation**: Help make docs clearer
- **Add examples**: Create tutorials and examples
- **Write tests**: Improve test coverage
- **Review code**: Comment on pull requests

## Getting Started

### 1. Fork and Clone

```bash
# Fork the repository on GitHub, then:
git clone https://github.com/YOUR_USERNAME/pypopart.git
cd pypopart
```

### 2. Set Up Development Environment

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode with dev dependencies
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

### 3. Create a Branch

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/bug-description
```

## Development Workflow

### Code Style

PyPopART uses:
- **ruff** for linting and formatting
- **mypy** for type checking
- **numpydoc** style docstrings

Format your code:
```bash
# Auto-format code
ruff format src/ tests/

# Check for issues
ruff check src/ tests/

# Fix auto-fixable issues
ruff check --fix src/ tests/
```

Type checking:
```bash
mypy src/
```

### Writing Tests

All new code should include tests. We use pytest:

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_sequence.py

# Run with coverage
pytest --cov=pypopart --cov-report=html
```

Test file structure:
```python
"""
Test module for XYZ functionality.
"""

import pytest
from pypopart.module import Class


class TestClass:
    """Test cases for Class."""
    
    def test_basic_functionality(self):
        """Test basic use case."""
        obj = Class()
        assert obj.method() == expected_result
    
    def test_error_handling(self):
        """Test error conditions."""
        with pytest.raises(ValueError):
            Class(invalid_param)
```

### Documentation

#### Docstrings

Use NumPy-style docstrings:

```python
def function(param1: int, param2: str) -> bool:
    """
    Brief description of function.

    Longer description with more details about what the function does.

    Parameters
    ----------
    param1 : int
        Description of param1.
    param2 : str
        Description of param2.

    Returns
    -------
    bool
        Description of return value.

    Raises
    ------
    ValueError
        When param1 is negative.

    Examples
    --------
    >>> function(42, "test")
    True
    """
    pass
```

#### Documentation Files

- Update relevant `.md` files in `docs/`
- Add examples to appropriate tutorial notebooks
- Update API reference if adding new classes/functions

### Commit Messages

Write clear commit messages:

```
Short summary (50 chars or less)

More detailed explanation if needed. Wrap at 72 characters.
Explain what changed and why, not how.

- Bullet points are okay
- Use present tense: "Add feature" not "Added feature"
- Reference issues: "Fixes #123"
```

Examples:
```
Add Tamura-Nei distance calculation

Implement the Tamura-Nei model for calculating evolutionary
distances, accounting for GC content and transition/transversion
rate differences.

Fixes #45
```

## Pull Request Process

### Before Submitting

1. **Tests pass**: `pytest`
2. **Code is formatted**: `ruff format src/ tests/`
3. **No linting errors**: `ruff check src/ tests/`
4. **Type checking passes**: `mypy src/`
5. **Documentation updated**: Relevant docs and docstrings
6. **CHANGELOG updated**: Add entry to `docs/changelog.md`

### Submitting

1. Push your branch to your fork
2. Open a pull request against `main`
3. Fill out the pull request template
4. Link relevant issues

### Pull Request Template

```markdown
## Description
Brief description of changes.

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Breaking change
- [ ] Documentation update

## Checklist
- [ ] Tests added/updated
- [ ] Documentation updated
- [ ] Code formatted with ruff
- [ ] Type hints added
- [ ] CHANGELOG updated

## Related Issues
Fixes #(issue number)
```

## Code Review Process

Maintainers will review your PR and may:
- Request changes
- Ask questions
- Suggest improvements

Please respond promptly to feedback. PRs inactive for 30 days may be closed.

## Reporting Bugs

### Before Reporting

1. Check existing issues
2. Try the latest version
3. Verify it's a PyPopART bug (not a dependency issue)

### Bug Report Template

```markdown
**Description**
Clear description of the bug.

**To Reproduce**
Steps to reproduce:
1. Load file '...'
2. Run command '...'
3. See error

**Expected Behavior**
What should have happened.

**Actual Behavior**
What actually happened.

**Environment**
- PyPopART version: 
- Python version:
- OS:

**Additional Context**
Error messages, screenshots, etc.
```

## Suggesting Features

Feature requests are welcome! Please provide:

1. **Use case**: Why is this feature needed?
2. **Proposed solution**: How should it work?
3. **Alternatives considered**: Other approaches?
4. **Examples**: Similar features in other tools?

## Release Process

(For maintainers)

1. Update version in `pyproject.toml`
2. Update `CHANGELOG.md`
3. Create release on GitHub
4. Build and upload to PyPI

## Code of Conduct

### Our Standards

- Be respectful and inclusive
- Welcome newcomers
- Accept constructive criticism
- Focus on what's best for the community
- Show empathy

### Unacceptable Behavior

- Harassment or discrimination
- Trolling or insulting comments
- Personal or political attacks
- Publishing others' private information

### Enforcement

Violations may result in:
1. Warning
2. Temporary ban
3. Permanent ban

Report issues to: [project maintainer email]

## Questions?

- Ask on GitHub Discussions
- Open an issue with the "question" label
- Check the [FAQ](docs/faq.md)

## License

By contributing, you agree that your contributions will be licensed under the GNU General Public License v3.0 or later.

## Recognition

Contributors are recognized in:
- `AUTHORS.md` file
- Release notes
- Documentation

Thank you for contributing to PyPopART! 🌳
