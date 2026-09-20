import functools
import http.server
import os
import tempfile
import threading
import unittest
from pathlib import Path

from visualization.bundle import build
from visualization.fixture import make_fixture


@unittest.skipUnless(os.environ.get("VIZ_BROWSER_TEST") == "1", "Opt in with VIZ_BROWSER_TEST=1")
class BrowserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from playwright.sync_api import sync_playwright
        cls.tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.tmp.name)
        cls.source = make_fixture(cls.root / "fixture")
        cls.output = cls.root / "export"
        build(cls.source, cls.output)
        handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(cls.output))
        cls.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.url = f"http://127.0.0.1:{cls.server.server_port}"
        cls.playwright = sync_playwright().start()
        cls.browser = cls.playwright.chromium.launch(args=["--enable-unsafe-swiftshader"])

    @classmethod
    def tearDownClass(cls):
        cls.browser.close()
        cls.playwright.stop()
        cls.server.shutdown()
        cls.server.server_close()
        cls.tmp.cleanup()

    def setUp(self):
        self.context = self.browser.new_context(viewport={"width": 1440, "height": 1000}, reduced_motion="reduce")
        self.external = []
        self.errors = []
        self.context.route("**/*", self.route)
        self.page = self.context.new_page()
        self.page.on("pageerror", lambda error: self.errors.append(str(error)))

    def route(self, route):
        if route.request.url.startswith(self.url) or route.request.url.startswith(("data:", "file:")):
            route.continue_()
        else:
            self.external.append(route.request.url)
            route.abort()

    def tearDown(self):
        self.context.close()

    def open(self):
        self.page.goto(self.url)
        self.page.locator("#welcome").wait_for()

    def test_structure_modes_and_filters(self):
        self.open()
        self.page.locator('.example-link').click()
        self.page.locator('[data-tab="structure"]').click()
        self.assertEqual(self.page.locator(".viewer canvas").count(), 1)
        self.assertEqual(self.page.locator('[data-view]').count(), 2)
        self.page.locator('[data-view="binding_site"]').click()
        self.assertEqual(self.page.locator('[data-view="binding_site"]').get_attribute("aria-pressed"), "true")
        self.page.locator('[data-view="candidate"]').click()
        self.page.locator('[data-tab="candidates"]').click()
        self.page.locator('.table-button[data-candidate="fixture-1"]').first.click()
        self.assertIn("50%", self.page.locator("#inspector").inner_text())
        self.page.get_by_text("View residue coverage", exact=True).click()
        self.page.locator('.heat-cell[data-residue="toy:A:7"]').first.click()
        self.assertEqual(self.page.locator("#candidate-label").inner_text(), "Fixture beta")
        self.assertIn("Toy A:7", self.page.locator("#residue-selection").inner_text())
        self.page.locator('[data-tab="candidates"]').click()
        self.page.locator("#novelty").select_option("known")
        self.assertNotIn("Fixture alpha", self.page.locator("#chapter-content").inner_text())
        self.page.locator("#rejected").check()
        self.assertIn("Fixture gamma", self.page.locator("#chapter-content").inner_text())
        self.page.locator("#modality").select_option("biologic")
        self.assertIn("No candidates match", self.page.locator("#chapter-content").inner_text())
        self.assertEqual(self.external, [])
        self.assertEqual(self.errors, [])

    def test_chat_preview_and_reset(self):
        self.open()
        self.page.locator('[data-prompt="Colorectal cancer"]').click()
        self.assertEqual(self.page.locator("#prompt").input_value(), "Colorectal cancer")
        self.page.locator("#prompt").press("Enter")
        self.assertIn("isn’t connected", self.page.locator(".assistant-message").inner_text())
        self.assertFalse(self.page.locator("#result").is_visible())
        self.page.locator(".inline-action").click()
        for tab in ("candidates", "structure", "sources"):
            self.page.locator(f'[data-tab="{tab}"]').click()
            self.assertIn("Synthetic example", self.page.locator("#origin-banner").inner_text())
        self.page.locator("#new-chat").click()
        self.assertTrue(self.page.locator("#welcome").is_visible())
        self.assertFalse(self.page.locator("#result").is_visible())
        self.assertEqual(self.page.locator("#conversation").inner_text(), "")
        self.assertEqual(self.external, [])
        self.assertEqual(self.errors, [])

    def test_no_webgl_fallback(self):
        self.context.add_init_script("HTMLCanvasElement.prototype.getContext = function() { return null; };")
        self.open()
        self.page.locator('.example-link').click()
        self.page.locator('[data-tab="structure"]').click()
        self.assertIn("3D unavailable", self.page.locator("#viewer-status").inner_text())
        self.page.locator('[data-tab="candidates"]').click()
        self.assertIn("Fixture alpha", self.page.locator("#chapter-content").inner_text())
        self.assertEqual(self.errors, [])

    def test_report_file_offline(self):
        self.context.set_offline(True)
        self.page.goto((self.output / "report.html").as_uri())
        self.assertIn("SYNTHETIC FIXTURE", self.page.locator("body").inner_text())
        self.assertIn("not evaluated", self.page.locator("body").inner_text().lower())
        self.assertEqual(self.external, [])

    def test_mobile_no_document_overflow(self):
        self.page.set_viewport_size({"width": 390, "height": 844})
        self.open()
        self.page.locator('.example-link').click()
        overflow = self.page.evaluate("document.documentElement.scrollWidth > window.innerWidth")
        self.assertFalse(overflow)
        self.assertEqual(self.errors, [])


if __name__ == "__main__":
    unittest.main()
