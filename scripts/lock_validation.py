"""Validate locked package identities and required versions, not a stale count."""
from collections.abc import Mapping


def locked_packages_valid(packages, required_versions):
    if not isinstance(packages, Mapping) or not packages:
        return False
    if any(
        not isinstance(record, Mapping)
        or record.get("Package") != name
        or not isinstance(record.get("Version"), str)
        or not record["Version"].strip()
        for name, record in packages.items()
    ):
        return False
    return all(
        packages.get(name, {}).get("Version") == version
        for name, version in required_versions.items()
    )
