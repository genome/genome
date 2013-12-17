-- Deploy web_nomenclature_permission
-- requires: web_nomenclature

BEGIN;

REVOKE ALL ON TABLE web.nomenclature FROM PUBLIC;
REVOKE ALL ON TABLE web.nomenclature FROM genome;
GRANT ALL ON TABLE web.nomenclature TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE web.nomenclature TO "gms-user";

COMMIT;
