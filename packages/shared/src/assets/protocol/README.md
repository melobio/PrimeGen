# Schema Bump Checklist

1. Create new JSON schema in `shared-data/protocol/schemas/${schemaVersion}.json`
2. Create TS types for new schema in `shared-data/protocol/types`
3. Create or modify Python types as necessary in both `python/opentrons_shared_data/protocol/dev_types` and in the python api `dev_types` in `api/src/opentrons/protocol_api/dev_types.py`
4. Create new executor in api to handle new schema. Be sure to include unit tests
5. Add new test in `api/tests/opentrons/test_execute.py` for the schema
6. Add new test in `shared-data/js/__tests__/protocolSchemaV${schemaVersion}.test.js` to ensure new schema definition functions as intended
7. Add new schema fixture to `shared-data/protocol/fixtures/${schemaVersion}/${schemaVersion}.json` that the above test will read
8. Update `MAX_SUPPORTED_VERSION` in `api/src/opentrons/protocols/parse.py`

**IMPORTANT**: Be sure to clearly communicate the schema bump to the software + QA teams. This is a strict contract that should be honored by both the frontend and robot stack.
