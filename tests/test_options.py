import unittest
from unittest.mock import Mock, mock_open, patch

from ppfit.options import read_options

class OptionsTest( unittest.TestCase ):

    def test_read_option( self ):
        example_file = """key: value"""
        with patch( 'builtins.open', mock_open( read_data=example_file ), create=True ) as m:
            self.assertEqual( read_options( 'filename' ), { 'key': 'value' } )


if __name__ == '__main__':
    unittest.main()
