import re
import yaml
from rocrate.rocrate import ROCrate


def extract_doi_parts(doi_string):
    """
    Extract and clean a DOI from arbitrary text (plain DOI string, full URL, etc.).

    Parameters:
    - doi_string (str): Any string that may contain a DOI.

    Returns:
    The cleaned DOI string (e.g. '10.1234/example'), or
    'No valid DOI found in the input string.' if no DOI is found.
    """
    doi_pattern = re.compile(r'(10\.[0-9]+/[^ \s]+)')
    match = doi_pattern.search(doi_string)

    if match:
        doi = match.group(1)
        # Strip trailing punctuation that may have been accidentally included
        doi = re.sub(r'[\s,.:;|\/\?:@&=+\$,]+$', '', doi)
        return doi
    else:
        return "No valid DOI found in the input string."


def format_citation(crate):
    """
    Build an APA-style citation string from an ROCrate object.

    Parameters:
    - crate (ROCrate): A loaded ROCrate instance.

    Returns:
    A citation string in the form:
      'Surname, I., & Surname, I. (Year). Title [Data set]. Publisher.
       https://doi.org/suffix'
    or an error string if the root entity is missing.
    """
    root = crate.root_dataset
    if root is None:
        return "Error: Root data entity not found."

    title = root.get("name") or "No title available"

    # identifier may be a string, a list, or absent
    identifier = root.get("identifier")
    if isinstance(identifier, list):
        doi = identifier[0] if identifier and identifier[0] else "No DOI available"
    elif isinstance(identifier, str) and identifier:
        doi = identifier
    else:
        doi = "No DOI available"

    # Extract year from datePublished
    date_published = (root.get("datePublished") or "")[:4]

    # Resolve publisher names
    publisher_refs = root.get("publisher") or []
    if not isinstance(publisher_refs, list):
        publisher_refs = [publisher_refs]
    publisher_names = []
    for ref in publisher_refs:
        entity = crate.dereference(ref["@id"]) if isinstance(ref, dict) else crate.dereference(ref.id)
        if entity:
            publisher_names.append(entity.get("name") or "No publisher available")
    publisher_str = ", ".join(publisher_names) if publisher_names else "No publisher available"

    # Resolve and format author names
    creator_refs = root.get("creator") or []
    if isinstance(creator_refs, dict):
        creator_refs = [creator_refs]
    author_names = []
    for ref in creator_refs:
        entity = crate.dereference(ref["@id"]) if isinstance(ref, dict) else crate.dereference(ref.id)
        if entity:
            surname = entity.get("familyName") or ""
            given = entity.get("givenName") or ""
            initial = given[0] if given else ""
            author_names.append(f"{surname}, {initial}.")

    if len(author_names) > 1:
        authors_str = ", ".join(author_names[:-1]) + f", & {author_names[-1]}"
    else:
        authors_str = "".join(author_names)

    # Use just the suffix portion of the DOI for the URL
    doi_suffix = doi.split("/")[-1] if doi != "No DOI available" else doi
    citation = (
        f"{authors_str} ({date_published}). {title} [Data set]. "
        f"{publisher_str}. https://doi.org/{doi_suffix}"
    )
    return citation


def ro_crate_to_cff(crate):
    """
    Convert an ROCrate object to a CFF-formatted YAML string.

    Parameters:
    - crate (ROCrate): A loaded ROCrate instance.

    Returns:
    A YAML string ready to write directly to CITATION.cff, or an error string
    if the root entity is missing.
    """
    root = crate.root_dataset
    if root is None:
        return "Error: Root data entity not found."

    title = root.get("name") or "No title available"
    version = root.get("version") or "1.0"

    identifier = root.get("identifier")
    if isinstance(identifier, list):
        doi = identifier[0] if identifier else "No DOI available"
    else:
        doi = identifier or "No DOI available"

    date_released = (root.get("datePublished") or "").split("T")[0]
    url = root.get("url") or "No URL provided"

    # Resolve authors
    creator_refs = root.get("creator") or []
    if isinstance(creator_refs, dict):
        creator_refs = [creator_refs]

    author_list = []
    for ref in creator_refs:
        if isinstance(ref, dict):
            author_id = ref.get("@id")
        else:
            author_id = ref.id

        if author_id is None:
            print(f"No '@id' found for author reference: {ref}")
            continue

        entity = crate.dereference(author_id)
        if entity:
            author_list.append({
                "family-names": entity.get("familyName") or "",
                "given-names": entity.get("givenName") or "",
                "orcid": author_id,
            })
        else:
            print(f"Could not dereference author entity: {author_id}")

    cff_dict = {
        "cff-version": "1.2.0",
        "message": "If you use this model, please cite it as below.",
        "authors": author_list,
        "title": title,
        "version": version,
        "doi": doi,
        "date-released": date_released,
        "url": url,
        "type": "dataset",
    }

    return yaml.dump(cff_dict, sort_keys=False, default_flow_style=False)
