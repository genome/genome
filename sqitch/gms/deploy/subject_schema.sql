-- Deploy subject_schema
-- requires: empty_db

BEGIN;

CREATE SCHEMA subject AUTHORIZATION genome;

COMMIT;
