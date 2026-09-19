import unittest
from lock_validation import locked_packages_valid


class LockValidationTests(unittest.TestCase):
    def test_additional_packages_do_not_break_the_chapter_contract(self):
        packages = {"A": {"Package": "A", "Version": "1.0"}}
        self.assertTrue(locked_packages_valid(packages, {"A": "1.0"}))
        packages["B"] = {"Package": "B", "Version": "2.0"}
        self.assertTrue(locked_packages_valid(packages, {"A": "1.0"}))

    def test_missing_wrong_version_and_malformed_records_fail(self):
        self.assertFalse(locked_packages_valid({}, {"A": "1.0"}))
        self.assertFalse(locked_packages_valid({"B": {"Package": "B", "Version": "1.0"}}, {"A": "1.0"}))
        self.assertFalse(locked_packages_valid({"A": {"Package": "A", "Version": "2.0"}}, {"A": "1.0"}))
        self.assertFalse(locked_packages_valid({"A": {"Package": "B", "Version": "1.0"}}, {"A": "1.0"}))


if __name__ == '__main__':
    unittest.main()
