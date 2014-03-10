-- Deploy config_analysismenu_item
-- requires: config_schema

BEGIN;
CREATE TABLE IF NOT EXISTS config.analysismenu_item (
  id character varying(64) PRIMARY KEY,
  created_at timestamp(6) without time zone NOT NULL,
  updated_at timestamp(6) without time zone NOT NULL,
  created_by character varying(255) NOT NULL,
  name text NOT NULL,
  file_path text NOT NULL
);

REVOKE ALL ON TABLE config.analysismenu_item FROM PUBLIC;
GRANT ALL ON TABLE config.analysismenu_item TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.analysismenu_item TO "gms-user";
COMMIT;
