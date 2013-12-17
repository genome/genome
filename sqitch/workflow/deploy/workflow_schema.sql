-- Deploy workflow_schema
-- requires: empty_db

BEGIN;

CREATE SCHEMA workflow AUTHORIZATION genome;

COMMIT;
