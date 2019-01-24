# PopBio-interchange-format
Simple Abundance Format (SAF) for VectorBase PopBio incidence/abundance datasets from arthropod surveillance studies. 20 column format to capture the minimum necessary metadata for a collection event. 



Description of Simple Abundance Format (SAF)
--------------------------------------------
Field Name |Format|Requirement|Details
-----------|------|-------|-----------
collection_ID|string|Mandatory|Identifier for collection event e.g. ABC_2018_collection_00001
sample_ID|string|Mandatory|Identifier for sample event e.g. ABC_2018_sample_00001
collection_start_date|ISO 8601 date format (YYYY-MM-DD)|Mandatory|Date at which traps deployed
collection_end_date|ISO 8601 date format (YYYY-MM-DD)|Mandatory|Date at which traps collected and specimens removed for processing
trap_ID|string|Advisory|Internal trap ID for the collection, may be used as part of processing for distinct collection events
GPS_latitude|GPS decimal degrees|Mandatory|Latitude for collection site, max. 6 decimal places
GPS_longitude|GPS decimal degrees|Mandatory|Longitude for collection site, max. 6 decimal places
GPS_qualifier|string|Optional|controlled vocabulary/ontology term to describe source and precision of GPS coordinates
location_description|string|Optional|Collection location description e.g. Orlando
location_ADM2|string|Optional|Administrative level 2 for collection e.g. Orange County
location_ADM1|string|Optional|Administrative level 1 for collection e.g. Florida
location_country|string|Optional|Country in which collection occurred e.g. United States of America
trap_type|string|Mandatory|Trap type e.g. CDC light trap, New Jersey Trap
attractant|string|Advisory|List of attractants used in the trap e.g. CO2, light
trap_number|integer|Mandatory|Number of traps deployed (Default is 1)
trap_duration|integer|Mandatory|Number of nights/days trap was deployed (Default is 1)
species|string|Species|binomial species name for collected specimens
species_identification_method|string|Mandatory|Protocol for asserting species identification
developmental_stage|string|Mandatory|developmental stage e.g. adult
sex|string|Mandatory|sex of specimens e.g. female/male/mixed
sample_count|integer|Mandatory|count of specimens from collection
collection_comment|string|Optional|free text comment about the collection site or event
sample_comment|string|Optional|free text comment about the sample material
species_comment|string|Optional|free text comments about the species identification process
