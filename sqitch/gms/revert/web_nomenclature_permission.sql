-- Revert web_nomenclature_permission

BEGIN;

REVOKE ALL ON TABLE web.nomenclature FROM "gms-user";

COMMIT;
