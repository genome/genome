-- Deploy instrument_fragment_library_permission
-- requires: instrument_fragment_library

BEGIN;

REVOKE ALL ON TABLE instrument.fragment_library FROM PUBLIC;
REVOKE ALL ON TABLE instrument.fragment_library FROM genome;
GRANT ALL ON TABLE instrument.fragment_library TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE instrument.fragment_library TO "gms-user";

COMMIT;
