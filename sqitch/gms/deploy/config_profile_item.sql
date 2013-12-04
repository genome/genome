-- Deploy config_profile_item
-- requires: config_analysis_project
-- requires: config_analysismenu_item

BEGIN;
CREATE TABLE IF NOT EXISTS config.profile_item (
  id character varying(64) PRIMARY KEY,
  created_at timestamp(6) without time zone NOT NULL,
  updated_at timestamp(6) without time zone NOT NULL,
  created_by character varying(255) NOT NULL,
  analysismenu_item_id character varying(64) REFERENCES config.analysismenu_item(id),
  analysis_project_id character varying(64) NOT NULL REFERENCES config.analysis_project(id)
);

REVOKE ALL ON TABLE config.profile_item FROM PUBLIC;
GRANT ALL ON TABLE config.profile_item TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.profile_item TO "gms-user";
COMMIT;
