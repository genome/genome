-- Deploy web_nomenclature_field_permission
-- requires: web_nomenclature_field

BEGIN;

REVOKE ALL ON TABLE web.nomenclature_field FROM PUBLIC;
REVOKE ALL ON TABLE web.nomenclature_field FROM genome;
GRANT ALL ON TABLE web.nomenclature_field TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE web.nomenclature_field TO "gms-user";

COMMIT;
