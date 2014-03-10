-- Deploy instrument_data_permission
-- requires: instrument_data

BEGIN;

REVOKE ALL ON TABLE instrument.data FROM PUBLIC;
REVOKE ALL ON TABLE instrument.data FROM genome;
GRANT ALL ON TABLE instrument.data TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE instrument.data TO "gms-user";

COMMIT;
