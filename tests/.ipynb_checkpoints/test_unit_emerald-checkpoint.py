from emeraldbgc import _cli
import pytest


def testCliHelp():
    with pytest.raises(SystemExit) as e:
        _cli.main(['--help'])
    assert e.value.code == 0

def testCliError():
    with pytest.raises(SystemExit) as e:
        _cli.main(['invalid'])
    assert e.value.code == 2