-- Deploy config_analysis_project_model_bridge
-- requires: config_analysis_project
-- requires: model_model

BEGIN;

CREATE TABLE IF NOT EXISTS config.analysis_project_model_bridge (
    id character varying(64) NOT NULL PRIMARY KEY,
    created_by character varying(255) NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    analysis_project_id character varying(64) NOT NULL REFERENCES config.analysis_project(id),
    model_id character varying(64) NOT NULL UNIQUE REFERENCES model.model(genome_model_id)
);

REVOKE ALL ON TABLE config.analysis_project_model_bridge FROM PUBLIC;
REVOKE ALL ON TABLE config.analysis_project_model_bridge FROM genome;
GRANT ALL ON TABLE config.analysis_project_model_bridge TO genome;
GRANT SELECT,INSERT,DELETE,UPDATE ON TABLE config.analysis_project_model_bridge TO "gms-user";

COMMIT;
