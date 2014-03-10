-- Revert subject_pairing

BEGIN;

DROP TABLE IF EXISTS subject.pairing;

COMMIT;
