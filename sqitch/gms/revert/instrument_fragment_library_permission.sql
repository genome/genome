-- Revert instrument_fragment_library_permission

BEGIN;

REVOKE ALL ON TABLE instrument.fragment_library FROM "gms-user";

COMMIT;
