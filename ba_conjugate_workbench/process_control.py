"""Process lifecycle helpers for the local Workbench server."""

from __future__ import annotations

import signal
import threading


def schedule_workbench_shutdown(delay_seconds: float = 1.5) -> threading.Timer:
    """Ask Streamlit to stop after the current UI response has been rendered.

    Streamlit installs a SIGTERM handler that shuts down its server, runtime, and
    worker threads. Raising the signal from a short-lived timer keeps this
    button scoped to the current Workbench process; a ``mamba run`` or
    ``conda run`` parent then exits normally when the server process finishes.
    """
    timer = threading.Timer(delay_seconds, signal.raise_signal, (signal.SIGTERM,))
    timer.daemon = True
    timer.start()
    return timer
