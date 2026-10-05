import signal

import ba_conjugate_workbench.process_control as process_control


def test_schedule_workbench_shutdown_uses_daemon_timer(monkeypatch):
    created = {}

    class FakeTimer:
        def __init__(self, interval, function, args):
            created["interval"] = interval
            created["function"] = function
            created["args"] = args
            self.daemon = False

        def start(self):
            created["started"] = True

    monkeypatch.setattr(process_control.threading, "Timer", FakeTimer)

    timer = process_control.schedule_workbench_shutdown(delay_seconds=2.0)

    assert created == {
        "interval": 2.0,
        "function": signal.raise_signal,
        "args": (signal.SIGTERM,),
        "started": True,
    }
    assert timer.daemon is True
