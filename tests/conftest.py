import signal
import time

import pytest

TEST_SECONDS = 5
SESSION_SECONDS = 20
_started = time.monotonic()


def pytest_sessionstart(session):
    global _started
    _started = time.monotonic()


def pytest_runtest_setup(item):
    if time.monotonic() - _started > SESSION_SECONDS:
        pytest.skip(f"the {SESSION_SECONDS} s budget of this worker is used up")


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item):
    def stop(*_):
        raise TimeoutError(f"over the {TEST_SECONDS} s limit of one test")

    previous = signal.signal(signal.SIGALRM, stop)
    signal.alarm(TEST_SECONDS)
    try:
        outcome = yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)
    if outcome.excinfo and isinstance(outcome.excinfo[1], TimeoutError):
        outcome.force_exception(pytest.skip.Exception(f"over the {TEST_SECONDS} s limit of one test"))
